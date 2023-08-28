#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
//#include <time.h>
long long time(int);
#include <complex.h>

#include "kspace_filter.h"
#include "noise.h"
#include "stdbool.h"


#include "calib/calib.h"
#include "calib/estvar.h"

#include "num/multind.h"
#include "num/flpmath.h"
#include "num/fft.h"
#include "num/init.h"
#include "num/ops_p.h"
#include "num/ops.h"
#include "num/iovec.h"

#include "iter/misc.h"
#include "iter/monitor.h"

#include "linops/linop.h"
#include "linops/fmac.h"
#include "linops/someops.h"

#include "noncart/nufft.h"

#include "sense/recon.h"
#include "sense/model.h"
#include "sense/optcom.h"
#include "sense/pocs.h"

#include "misc/debug.h"
#include "misc/memcfl.h"
#include "misc/mri.h"
#include "misc/utils.h"
#include "misc/mmio.h"
#include "misc/misc.h"
#include "misc/opts.h"

#include "grecon/optreg.h"
#include "grecon/italgo.h"


//#include "_kiss_fft_guts.h"
#include "params.h"

#define NA_SCALE 140

const float na_vol = 0.7;
const float na_t2fr = 60;
const float na_mm = 140;
const float na_t1 = 39.2;

const char* img_sim = "img_sim.mem";
const char* kspace_sim = "kspace_sim.mem";
const char* sen_sim = "sen_sim.mem";
const char* kspace_cs = "kspace_cs.mem";
const char* img_cs = "img_cs.mem";
const char* kspace_filt = "kspace_filt.mem";
const char* img_filt = "img_filt.mem";

const static int D = 16;

int main_ecalib(int argc, const char* argv[]);
int main_pics(int argc, const char* argv[]);
int main_fft(int argc, const char* argv[]);
int main_zeros(int argc, const char* argv[]);


static inline void normalizecimage(complex float* cimage, struct Params *p) {
    double max = 0;
    switch(p->sequence) {
        case Na:
        case NaTQ:
        case NaSQ:
            max = creal(cimage[0]);
            for(int i=0;i<p->ydim*p->zdim*p->xdim;i++) {
                if (max < creal(cimage[i])) {
                    max = creal(cimage[i]);
                }
            }
            for(int i=0;i<p->ydim*p->zdim*p->xdim;i++) {
                cimage[i]=creal(cimage[i])/max+0*I;
            }
            break;
        default:
            break;
    }
}

int call_count = 0;
static inline double simVoxel(struct Params *p, int pos, struct Dataset *ds) {
    double te, tr, ti, fa, tau1, tau2, te_end, te_step, e_tr_t1, e_tr_t2, cfa, sfa, s, m;
    double pd, t1, t2, t2s, t2f, ex_frac;
    call_count += 1;
    s = 0;
    if(pos >= ds->k_xdim*ds->k_ydim*ds->k_zdim) {
        fprintf(stdout, "%i\n", pos);
        return 0;
    }
    switch(p->sequence) {
        case SE:
            te = p->s_params[0];
            tr = p->s_params[1];
            pd = ds->pd[pos];
            t1 = ds->t1[pos];
            t2 = ds->t2[pos];
            s = fabs(pd*exp(-te/t2)*(1.0-exp(-tr/t1)));
            break;
        case IR:
            te = p->s_params[0];
            tr = p->s_params[1];
            ti = p->s_params[2];
            pd = ds->pd[pos];
            t1 = ds->t1[pos];
            t2 = ds->t2[pos];
            s = fabs(pd*(1.0-2.0*exp(-ti/t1)+exp(-tr/ti))*exp(-te/t2));
            break;
        case bSSFP:
            te = p->s_params[0];
            tr = p->s_params[1];
            fa = p->s_params[2] * M_PI / 180.0;
            pd = ds->pd[pos];
            t1 = ds->t1[pos];
            t2 = ds->t2[pos];
            e_tr_t1 = exp(-tr/t1);
            e_tr_t2 = exp(-tr/t2);
            sfa = sin(fa);
            cfa = cos(fa);
            s = fabs(pd*sfa*(1.0-e_tr_t1)/(1.0-(e_tr_t1-e_tr_t2)*cfa-e_tr_t1*e_tr_t2)*exp(-te/t2) );
            break;
        case pcbSSFP:
            te = p->s_params[0];
            tr = p->s_params[1];
            fa = p->s_params[2] * M_PI / 180.0;
            pd = ds->pd[pos];
            t1 = ds->t1[pos];
            t2 = ds->t2[pos];
            e_tr_t1 = exp(-tr/t1);
            e_tr_t2 = exp(-tr/t2);
            sfa = sin(fa);
            cfa = cos(fa);
            // https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fmrm.29302&file=mrm29302-sup-0001-supinfo.pdf
            break;
        case FISP:
            te = p->s_params[0];
            tr = p->s_params[1];
            fa = p->s_params[2] * M_PI / 180.0;
            pd = ds->pd[pos];
            t1 = ds->t1[pos];
            t2 = ds->t2[pos];
            t2s = ds->t2s[pos];
            e_tr_t1 = exp(-tr/t1);
            e_tr_t2 = exp(-tr/t2);
            sfa = sin(fa);
            cfa = cos(fa);
            s = fabs(pd*sfa/(1+cfa) * (1 - (e_tr_t1-cfa)*sqrt((1-e_tr_t2*e_tr_t2)/( (1-e_tr_t1*cfa)*(1-e_tr_t1*cfa)-e_tr_t2*e_tr_t2*(e_tr_t1-cfa)*(e_tr_t1-cfa) ) ) ) * exp(-te/t2s));
            break;
        case PSIF:
            te = p->s_params[0];
            tr = p->s_params[1];
            fa = p->s_params[2] * M_PI / 180.0;
            pd = ds->pd[pos];
            t1 = ds->t1[pos];
            t2 = ds->t2[pos];
            e_tr_t1 = exp(-tr/t1);
            e_tr_t2 = exp(-tr/t2);
            sfa = sin(fa);
            cfa = cos(fa);
            s = fabs(pd * sfa/(1+cfa)*(1-(1-e_tr_t1*cfa)*sqrt((1-e_tr_t2*e_tr_t2) /( (1-e_tr_t1*cfa)*(1-e_tr_t1*cfa)-e_tr_t2*e_tr_t2*(e_tr_t1-cfa)*(e_tr_t1-cfa) ))) * exp(-te/t2));
            break;
        case SGRE:
            te = p->s_params[0];
            tr = p->s_params[1];
            fa = p->s_params[2] * M_PI / 180.0;
            pd = ds->pd[pos];
            t1 = ds->t1[pos];
            t2s = ds->t2s[pos];
            sfa = sin(fa);
            cfa = cos(fa);
            e_tr_t1 = exp(-tr/t1);
            s = fabs(pd * (1-e_tr_t1)*sfa/(1-e_tr_t1*cfa) * exp(-te/t2s));
            break;
        case Na:
            te = p->s_params[0];
            tr = p->s_params[1];
            ex_frac = ds->na_ex_frac[pos];
            t1 = ds->na_t1[pos];
            t2s = ds->na_t2s[pos];
            t2f = ds->na_t2f[pos];
            pd = ds->na_mm[pos];
            if(t2f==0) {
                s = fabs((na_vol-ex_frac)*pd* (1-exp(-tr/t1)) * (0.6 + 0.4*exp(-te/t2s)) + ex_frac*na_mm* (1-exp(-tr/na_t1))*exp(-te/na_t2fr)) / NA_SCALE;
            } else {
            s = fabs((na_vol-ex_frac)*pd* (1-exp(-tr/t1)) * (0.6*exp(-te/t2f) + 0.4*exp(-te/t2s)) + ex_frac*na_mm* (1-exp(-tr/na_t1))*exp(-te/na_t2fr)) / NA_SCALE;
            }
            break;
        case NaSQ:
            fa = p->s_params[0] * M_PI / 180;
            tau1 = p->s_params[1];
            te = p->s_params[2];
            te_end = p->s_params[3];
            te_step = p->s_params[4];
            t2s = ds->na_t2s[pos];
            t2f = ds->na_t2f[pos];
            pd = ds->na_mm[pos];
            sfa = sin(fa);
            s = 0;
            m = 0;
            for(;te<=te_end;te+=te_step) {
                m++;
                if(t2f == 0) {
                    s += fabs(pd*(exp(-(te+tau1)/t2s) + 1)*sfa);
                } else {
                    s += fabs(pd*(exp(-(te+tau1)/t2s) + exp(-(te+tau1)/t2f))*sfa);
                }
            }
            s /= m;
            break;
        case NaTQ:
            fa = p->s_params[0] * M_PI / 180;
            tau1 = p->s_params[1];
            tau2 = p->s_params[2];
            te = p->s_params[3];
            te_end = p->s_params[4];
            te_step = p->s_params[5];
            t2s = ds->na_t2s[pos];
            t2f = ds->na_t2f[pos];
            pd = ds->na_mm[pos];
            sfa = sin(fa);
            s = 0;
            m = 0;
            for(;te<=te_end;te+=te_step) {
                m++;
                if(t2f==0) {
                    s += fabs(pd*( (exp(-te/t2s) - 1) * (exp(-tau1/t2s)-1) * exp(-tau2/t2s) * (sfa*sfa*sfa*sfa*sfa) ));
                } else {
                    s += fabs(pd*( (exp(-te/t2s) - exp(-te/t2f)) * (exp(-tau1/t2s)-exp(-tau1/t2f)) * exp(-tau2/t2s) * (sfa*sfa*sfa*sfa*sfa) ));
                }
            }
            s /= m;
            break;
        case NaTQSQ:
            fa = p->s_params[0] * M_PI / 180;
            tau1 = p->s_params[1];
            tau2 = p->s_params[2];
            te = p->s_params[3];
            te_end = p->s_params[4];
            te_step = p->s_params[5];
            t2s = ds->na_t2s[pos];
            t2f = ds->na_t2f[pos];
            pd = ds->na_mm[pos];
            sfa = sin(fa);
            e_tr_t1 = 0;
            m = 0;
            for(;te<=te_end;te+=te_step) {
                m++;
                if(t2f==0) {
                    e_tr_t1 += fabs(pd*( (exp(-te/t2s) - 1) * (exp(-tau1/t2s)-1) * exp(-tau2/t2s) * (sfa*sfa*sfa*sfa*sfa) ));
                } else {
                    e_tr_t1 += fabs(pd*( (exp(-te/t2s) - exp(-te/t2f)) * (exp(-tau1/t2s)-exp(-tau1/t2f)) * exp(-tau2/t2s) * (sfa*sfa*sfa*sfa*sfa) ));
                }
            }
            e_tr_t1 /= m*NA_SCALE;
            fa = p->s_params[6] * M_PI / 180;
            tau1 = p->s_params[7];
            te = p->s_params[8];
            te_end = p->s_params[9];
            te_step = p->s_params[10];
            sfa = sin(fa);
            e_tr_t2 = 0;
            m = 0;
            for(;te<=te_end;te+=te_step) {
                m++;
                if(t2f == 0) {
                    e_tr_t2 += fabs(pd*(exp(-(te+tau1)/t2s) + 1)*sfa);
                } else {
                    e_tr_t2 += fabs(pd*(exp(-(te+tau1)/t2s) + exp(-(te+tau1)/t2f))*sfa);
                }
            }
            e_tr_t2 /= m*NA_SCALE;
            s = e_tr_t1 / e_tr_t2;
            break;
        case NaTQF:
            te = p->s_params[0];
            fa = p->s_params[1] * M_PI / 180;
            tau1 = p->s_params[2];
            tau2 = p->s_params[3];
            t2s = ds->na_t2s[pos];
            t2f = ds->na_t2f[pos];
            pd = ds->na_mm[pos];
            sfa = sin(fa);
            cfa = cos(fa);
            if(t2f==0) {
                s = fabs(pd*(exp(-te/t2s)-1) * (exp(-tau1/t2s)-1) * exp(-tau2/t2s) * (sfa*sfa*sfa*sfa*sfa));
            } else {
                s = fabs(pd*(exp(-te/t2s)-exp(-te/t2f)) * (exp(-tau1/t2s)-exp(-tau1/t2f)) * exp(-tau2/t2s) * (sfa*sfa*sfa*sfa*sfa));
            }
            break;
    }
    return s;
}

double ssim(complex float* img1, complex float* img2, int size) {
    double mu1, mu2, var1, var2, sig_co, c1, c2, L, k1, k2;
    k1 = 0.01;
    k2 = 0.03;
    //L = pow(2.0, sizeof(kiss_fft_scalar)*8)-1;
    L = 0;
    for(int i=0;i<size;i++) {
        if(creal(img1[i]) > L) {
            L = creal(img1[i]);
        }
    }
    c1 = (L*k1)*(L*k1);
    c2 = (L*k2)*(L*k2);
    mu1 = 0.0;
    mu2 = 0.0;
    for(int i=0;i<size;i++) {
        mu1 += creal(img1[i]);
        mu2 += creal(img2[i]);
    }
    mu1 /= size;
    mu2 /= size;
    var1 = 0.0;
    var2 = 0.0;
    sig_co = 0.0;
    for(int i=0;i<size;i++) {
        var1 += (mu1-creal(img1[i]))*(mu1-creal(img1[i]));
        var2 += (mu2-creal(img2[i]))*(mu2-creal(img2[i]));
        sig_co += (mu2-creal(img2[i]))*(mu1-creal(img1[i]));
    }
    var1 /= size;
    var2 /= size;
    sig_co /= size;
    double ssim = (2.0*mu1*mu2+c1)*(2.0*sig_co+c2)/((mu1*mu1+mu2*mu2+c1)*(var1+var2+c2));
    //fprintf(stdout, "ssim: %f\n", ssim);
    return ssim;
}

void debug_print_files() {
    struct memcfl* mem = get_list();
    debug_printf(DP_INFO, "Files:\n");
	while (NULL != mem) {
        //debug_printf(DP_INFO, "%i %i %i %i\n", mem->name, mem->refcount, mem->managed, mem->next);
        debug_printf(DP_INFO, "%s\n", mem->name);
		mem = mem->next;
	}
}

void remove_all_files() {
    struct memcfl* mem = get_list();
	while (NULL != mem) {
        mem->refcount = 0;
        //debug_printf(DP_INFO, "%s %i\n", mem->name, mem->refcount);
        memcfl_unlink(mem->name);
        mem = get_list();
	}
}

void sim_multi_coil(const char* in_file, const char* out_file, long cal_dims[D], long ncoils) {
    int N = D;
	long ksp_dims[N];
	complex float* in_data = memcfl_load(in_file, N, ksp_dims);

    md_copy_dims(DIMS, cal_dims, ksp_dims);
    cal_dims[COIL_DIM] = ncoils;
    complex float* coil_data = memcfl_create(out_file, D, cal_dims);

    for(int x=0;x<cal_dims[0];x++) {
        for(int y=0;y<cal_dims[1];y++) {
            for(int z=0;z<cal_dims[2];z++) {
                for(int c=0;c<cal_dims[COIL_DIM];c++) {
                    int opos = x + y*cal_dims[0] + z*cal_dims[0]*cal_dims[1] + c*cal_dims[0]*cal_dims[1]*cal_dims[2];
                    int ipos = x + y*cal_dims[0] + z*cal_dims[0]*cal_dims[1];
                    float w = 1;
                    switch(c) {
                        case 0:
                            w = (double)x/(double)cal_dims[0];
                            break;
                        case 1:
                            w = 1.0-(double)x/(double)cal_dims[0];
                            break;
                        case 2:
                            w = (double)y/(double)cal_dims[1];
                            break;
                        case 3:
                            w = 1.0-(double)y/(double)cal_dims[1];
                            break;
                    }
                    coil_data[opos] = in_data[ipos]*w;
                }
            }
        }
    }
    unmap_cfl(D, cal_dims, coil_data);
    unmap_cfl(D, ksp_dims, in_data);
}

void simulate(struct Params *p, struct Dataset *ds) {
    
    double xs,ys,zs, nx,ny,nz, d,s;
    int x,y,z, ipos,ds_pos;
    int xi,x_start,x_end,yi,y_start,y_end,zi,z_start,z_end;

    complex float *sim_img, *sim_kspace, *cs_kspace, *filt_img, *filt_kspace;

    char *flags = NULL;
    if(p->use_fft3) {
        flags = "7";
    } else {
        flags = "3";
    }

    long dims[D] = {p->xdim, p->ydim, p->zdim, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    remove_all_files();

    char sxdim[8];
    char sydim[8];
    char szdim[8];
    sprintf(sxdim, "%i", p->xdim);
    sprintf(sydim, "%i", p->ydim);
    sprintf(szdim, "%i", p->zdim); 
    const char* argv_zeros[] = {"zeros", "16", sxdim, sydim, szdim, "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "img_sim.mem"};
    main_zeros(19, argv_zeros);
    sim_img = memcfl_load(img_sim, D, dims);

    srand(time(0));

    call_count = 0;
    xs = (double)ds->k_xdim/(double)p->xdim/2.0;
    ys = (double)ds->k_ydim/(double)p->ydim/2.0;
    zs = (double)ds->k_zdim/(double)p->zdim/2.0;

    for(int i=0;i<p->zdim*p->ydim*p->xdim;i++) {
        sim_img[i] =  0+ 0*I;
    }

    p->ncoils = 1;
    for(z = 0; z<p->zdim; z++) {
        nz = (double)z*(double)ds->k_zdim/(double)p->zdim + zs;
        for(y = 0;y<p->ydim; y++) {
            ny = (double)y*(double)ds->k_ydim/(double)p->ydim + ys;
            for(x = 0;x<p->xdim; x++) {
                ipos = z*p->xdim*p->ydim + y*p->xdim + x;
                d = 0; s = 0;
                nx = (double)x*(double)ds->k_xdim/(double)p->xdim + xs;
                switch (p->nearest) {
                    case Nearest:
                        ds_pos = round(nz-0.5)*ds->k_xdim*ds->k_ydim+round(ny-0.5)*ds->k_xdim+round(nx-0.5);
                        s = simVoxel(p, ds_pos, ds);
                        break;
                    case KNearest:
                        ds_pos = floor(nz)*ds->k_xdim*ds->k_ydim+floor(ny)*ds->k_xdim+floor(nx);
                        s = simVoxel(p, ds_pos, ds);
                        ds_pos = floor(nz)*ds->k_xdim*ds->k_ydim+floor(ny)*ds->k_xdim+ceil(nx);
                        s += simVoxel(p, ds_pos, ds);
                        ds_pos = floor(nz)*ds->k_xdim*ds->k_ydim+ceil(ny)*ds->k_xdim+floor(nx);
                        s += simVoxel(p, ds_pos, ds);
                        ds_pos = floor(nz)*ds->k_xdim*ds->k_ydim+ceil(ny)*ds->k_xdim+ceil(nx);
                        s += simVoxel(p, ds_pos, ds);
                        ds_pos = ceil(nz)*ds->k_xdim*ds->k_ydim+floor(ny)*ds->k_xdim+floor(nx);
                        s += simVoxel(p, ds_pos, ds);
                        ds_pos = ceil(nz)*ds->k_xdim*ds->k_ydim+floor(ny)*ds->k_xdim+ceil(nx);
                        s += simVoxel(p, ds_pos, ds);
                        ds_pos = ceil(nz)*ds->k_xdim*ds->k_ydim+ceil(ny)*ds->k_xdim+floor(nx);
                        s += simVoxel(p, ds_pos, ds);
                        ds_pos = ceil(nz)*ds->k_xdim*ds->k_ydim+ceil(ny)*ds->k_xdim+ceil(nx);
                        s += simVoxel(p, ds_pos, ds);
                        s /= 8.0;
                        break;
                    case Average:
                    default:
                        x_start=floor(nx-xs);
                        x_end=ceil(x_start+2*xs);
                        y_start=floor(ny-ys);
                        y_end=ceil(y_start+2*ys);
                        z_start=floor(nz-zs);
                        z_end=ceil(z_start+2*zs);
                        for(xi=x_start;xi<x_end;xi+=1) {
                            for(yi=y_start;yi<y_end;yi+=1) {
                                for(zi=z_start;zi<z_end;zi+=1) {
                                    if(zi>=0&&zi<ds->k_zdim&&yi>=0&&yi<ds->k_ydim&&xi>=0&&xi<ds->k_xdim) {
                                        ds_pos = zi*ds->k_xdim*ds->k_ydim+yi*ds->k_xdim+xi;
                                        s += simVoxel(p, ds_pos, ds);
                                        d += 1.0;
                                    }
                                }
                            }
                        }
                        s /= d;
                        break;
                }
                sim_img[ipos] =  s+ 0*I;
            }
        }
    }

    normalizecimage(sim_img, p);

    addImageNoise(sim_img, p);

    memcfl_unmap(sim_img);

    const char* argv_fft[] = {"fft", "-u", flags, "img_sim.mem", "kspace_sim.mem"};
    main_fft(5, argv_fft);

    long *idims = NULL;
    sim_kspace = memcfl_load("kspace_sim.mem", D, idims);
    if (addKSpaceNoise(sim_kspace, p)) {
        //fprintf(stdout, "ifft kspace_sim img_sim\n");
        memcfl_unlink("img_sim.mem");
        const char* argv_ifft[] = {"fft", "-i", "-u", flags, "kspace_sim.mem", "img_sim.mem"};
        main_fft(6, argv_ifft);
    }
    memcfl_unmap(sim_kspace);
    
    if(p->use_cs) {
        long nc_dims[D];
        sim_multi_coil(img_sim, kspace_cs, nc_dims, 4);
        p->ncoils = 4;
        const char* argv_fft2[] = {"fft", "-u", flags, "kspace_cs.mem", "kspace_filt.mem"};
        main_fft(5, argv_fft2);
        filt_kspace = memcfl_load(kspace_filt, D, nc_dims);
        apply_kspace_filter(filt_kspace, filt_kspace, p);
        const char* argv_ifft2[] = {"fft", "-i", "-u", flags, "kspace_filt.mem", "img_filt.mem"};
        main_fft(6, argv_ifft2);
        cs_kspace = memcfl_load(kspace_cs, D, nc_dims);
        for(int i=0;i<io_calc_size(D, nc_dims, 1);i++) {
            cs_kspace[i] = filt_kspace[i];
        }
        memcfl_unmap(filt_kspace);
        memcfl_unmap(cs_kspace);
        //fprintf(stdout, "ecalib\n");
        const char* argv_ecalib[] = {"ecalib","kspace_filt.mem","sen_sim.mem"};
        main_ecalib(3, argv_ecalib);
        const char* argv_pics[] = {"pics", "-l1", "-r0.001", "kspace_cs.mem", "sen_sim.mem", "img_cs.mem"};
        main_pics(6, argv_pics);
    } else {
        sim_kspace = memcfl_load(kspace_sim, D, idims);
        if(apply_kspace_filter(sim_kspace, sim_kspace, p)) {
            //fprintf(stdout, "ifft kspace_sim img_sim\n");
            memcfl_unlink("img_sim.mem");
            const char* argv_ifft[] = {"fft", "-i", "-u", flags, "kspace_sim.mem", "img_sim.mem"};
            main_fft(6, argv_ifft);
        }
        memcfl_unmap(sim_kspace);
    }

    //fprintf(stdout, "voxel sim call count: %d overlap: %f\n", call_count, (double)call_count/(256.0*256.0*256.0));
}
