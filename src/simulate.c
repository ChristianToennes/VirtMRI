#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include "fft.h"
#include "fft2.h"
#include "fft3.h"
#include "tinycs.h"
#include "stdbool.h"

#include "_kiss_fft_guts.h"

typedef enum Sequence {
    SE, IR, bSSFP, FISP, PSIF, SGRE, Na, NaSQ, NaTQ, NaTQSQ, NaTQF
} sequence;

typedef enum Nearest {
    Nearest, KNearest, Average
} nearest;

typedef struct {
    sequence sequence;
    int n_params;
    float* s_params;
    int xdim;
    int ydim;
    int zdim;
    nearest nearest;
    bool use_cs;
    bool use_fft3;
    cs_params* cs_params;
} params;

params* make_params(sequence sequence, int n_params, float* s_params, int xdim, int ydim, int zdim, nearest nearest, bool use_cs, bool use_fft3, cs_params* cs_params) {
    params* p = (params*)malloc(sizeof(params));
    p->sequence = sequence;
    p->n_params = n_params;
    p->s_params = s_params;
    p->xdim = xdim;
    p->ydim = ydim;
    p->zdim = zdim;
    p->nearest = nearest;
    p->use_cs = use_cs;
    p->use_fft3 = use_fft3;
    p->cs_params = cs_params;
    return p;
}

typedef struct {
    float* pd;
    float* t1;
    float* t2;
    float* t2s;
    float* na_mm;
    float* na_t1;
    float* na_ex_frac;
    float* na_t2s;
    float* na_t2f;
    int k_xdim;
    int k_ydim;
    int k_zdim;
} dataset;

dataset* make_dataset(int kxdim, int kydim, int kzdim, float* pd, float* t1, float* t2, float* t2s, float* na_mm, float* na_t1, float* na_ex_frac, float* na_t2s, float* na_t2f) {
    dataset* ds = (dataset*)malloc(sizeof(dataset));
    ds->k_xdim = kxdim;
    ds->k_ydim = kydim;
    ds->k_zdim = kzdim;
    ds->pd = pd;
    ds->t1 = t1;
    ds->t2 = t2;
    ds->t2s = t2s;
    ds->na_mm = na_mm;
    ds->na_t1 = na_t1;
    ds->na_ex_frac = na_ex_frac;
    ds->na_t2s = na_t2s;
    ds->na_t2f = na_t2f;
    return ds;
}

void free_dataset(dataset* ds) {
    free(ds->pd);
    free(ds->t1);
    free(ds->t2);
    free(ds->t2s);
    free(ds->na_mm);
    free(ds->na_t1);
    free(ds->na_ex_frac);
    free(ds->na_t2s);
    free(ds->na_t2f);
    free(ds);
}
const float na_vol = 0.7;
const float na_t2fr = 60;

void addGaussianNoise(kiss_fft_cpx* image) {
    
}

double simVoxel(params* p, int pos, dataset* ds) {
    double te, tr, ti, fa, tau1, tau2, te_end, te_step, e_tr_t1, e_tr_t2, cfa, sfa, s;
    s = 0;
    switch(p->sequence) {
        case SE:
            te = p->s_params[0];
            tr = p->s_params[1];
            s = fabs((double)ds->pd[pos]*exp(-te/(double)ds->t2[pos])*(1.0-exp(-tr/(double)ds->t1[pos])));
            break;
        case IR:
            te = p->s_params[0];
            tr = p->s_params[1];
            ti = p->s_params[2];
            s = fabs((double)ds->pd[pos]*(1.0-2.0*exp(-ti/(double)ds->t1[pos])+exp(-tr/ti))*exp(-te/(double)ds->t2[pos]));
            break;
        case bSSFP:
            te = p->s_params[0];
            tr = p->s_params[1];
            fa = p->s_params[2] * M_PI / 180.0;
            e_tr_t1 = exp(-tr/ds->t1[pos]);
            e_tr_t2 = exp(-tr/ds->t2[pos]);
            sfa = sin(fa);
            cfa = cos(fa);
            s = fabs(ds->pd[pos]*sfa*(1.0-exp(-tr/ds->t1[pos]))/(1.0-(e_tr_t1+e_tr_t2)*cfa-e_tr_t1*e_tr_t2)*exp(-te/ds->t2[pos]) );
            break;
        case FISP:
            te = p->s_params[0];
            tr = p->s_params[1];
            fa = p->s_params[2] * M_PI / 180.0;
            e_tr_t1 = exp(-tr/ds->t1[pos]);
            e_tr_t2 = exp(-tr/ds->t2[pos]);
            sfa = sin(fa);
            cfa = cos(fa);
            s = fabs(ds->pd[pos] * sfa/(1+cfa) * (1 - (e_tr_t1-cfa)*sqrt((1-e_tr_t2*e_tr_t2)/( (1-e_tr_t1*cfa)*(1-e_tr_t1*cfa)-e_tr_t2*e_tr_t2*(e_tr_t1-cfa)*(e_tr_t1-cfa) ) ) ) * exp(-te/ds->t2s[pos]));
            break;
        case PSIF:
            te = p->s_params[0];
            tr = p->s_params[1];
            fa = p->s_params[2] * M_PI / 180.0;
            e_tr_t1 = exp(-tr/ds->t1[pos]);
            e_tr_t2 = exp(-tr/ds->t2[pos]);
            sfa = sin(fa);
            cfa = cos(fa);
            s = fabs(ds->pd[pos] * sfa/(1+cfa)*(1-(1-e_tr_t1*cfa)*sqrt((1-e_tr_t2*e_tr_t2) /( (1-e_tr_t1*cfa)*(1-e_tr_t1*cfa)-e_tr_t2*e_tr_t2*(e_tr_t1-cfa)*(e_tr_t1-cfa) ))) * exp(-te/ds->t2[pos]));
            break;
        case SGRE:
            te = p->s_params[0];
            tr = p->s_params[1];
            fa = p->s_params[2] * M_PI / 180.0;
            sfa = sin(fa);
            cfa = cos(fa);
            e_tr_t1 = exp(-tr/ds->t1[pos]);
            s = fabs(ds->pd[pos] * (1-e_tr_t1)*sfa/(1-e_tr_t1*cfa) * exp(-te/ds->t2s[pos]));
            break;
        case Na:
            te = p->s_params[0];
            tr = p->s_params[1];
            s = fabs((na_vol-ds->na_ex_frac[pos])* ds->na_mm[pos] * (1-exp(-tr/ds->na_t1[pos])) * (0.6*exp(-te/ds->na_t2f[pos]) + 0.4*exp(-te/ds->na_t2s[pos])) + ds->na_ex_frac[pos]*ds->na_mm[pos] * (1-exp(-tr/ds->na_t1[pos]))*exp(-te/na_t2fr))/140;
            break;
        case NaSQ:
            fa = p->s_params[0] * M_PI / 180;
            tau1 = p->s_params[1];
            tau2 = p->s_params[2];
            te = p->s_params[3];
            te_end = p->s_params[4];
            te_step = p->s_params[5];
            sfa = sin(fa);
            s = 0;
            for(;te<te_end;te++) {
                s += fabs(ds->na_mm[pos]*(exp(-(te+tau1)/ds->na_t2s[pos]) * exp(-(te+tau1)/ds->na_t2f[pos]))*sfa);
            }
            s /= (te_step*140);
            break;
        case NaTQ:
            fa = p->s_params[0] * M_PI / 180;
            tau1 = p->s_params[1];
            tau2 = p->s_params[2];
            te = p->s_params[3];
            te_end = p->s_params[4];
            te_step = p->s_params[5];
            sfa = sin(fa);
            s = 0;
            for(;te<te_end;te++) {
                s += fabs(ds->na_mm[pos] * ( (exp(-te/ds->na_t2s[pos]) - exp(-te/ds->na_t2f[pos])) * (exp(-tau1/ds->na_t2s[pos])-exp(-tau1/ds->na_t2f[pos])) * exp(-tau2/ds->na_t2s[pos]) ));
            }
            s /= (te_step*140);
            break;
        case NaTQSQ:
            fa = p->s_params[0] * M_PI / 180;
            tau1 = p->s_params[1];
            tau2 = p->s_params[2];
            te = p->s_params[3];
            te_end = p->s_params[4];
            te_step = p->s_params[5];
            sfa = sin(fa);
            e_tr_t1 = 0;
            for(;te<te_end;te++) {
                e_tr_t1 += fabs(ds->na_mm[pos] * ( (exp(-te/ds->na_t2s[pos]) - exp(-te/ds->na_t2f[pos])) * (exp(-tau1/ds->na_t2s[pos])-exp(-tau1/ds->na_t2f[pos])) * exp(-tau2/ds->na_t2s[pos]) ));
            }
            e_tr_t1 /= (te_step*140);
            fa = p->s_params[6] * M_PI / 180;
            tau1 = p->s_params[7];
            tau2 = p->s_params[8];
            te = p->s_params[9];
            te_end = p->s_params[10];
            te_step = p->s_params[11];
            sfa = sin(fa);
            e_tr_t2 = 0;
            for(;te<te_end;te++) {
                e_tr_t2 += fabs(ds->na_mm[pos]*(exp(-(te+tau1)/ds->na_t2s[pos]) * exp(-(te+tau1)/ds->na_t2f[pos]))*sfa);
            }
            e_tr_t2 /= (te_step*140);
            s = e_tr_t1 / e_tr_t2;
            break;
        case NaTQF:
            te = p->s_params[0];
            fa = p->s_params[1] * M_PI / 180;
            tau1 = p->s_params[2];
            tau2 = p->s_params[3];
            sfa = sin(fa);
            cfa = cos(fa);
            s = fabs(ds->na_mm[pos] * (exp(-te/ds->na_t2s[pos])-exp(-te/ds->na_t2f[pos])) * (exp(-tau1/ds->na_t2s[pos])-exp(-tau1/ds->na_t2f[pos])) * exp(-tau2/ds->na_t2s[pos]) * (sfa*sfa*sfa*sfa*sfa));
            break;
    }
    return s;
}

void simulate(params* p, kiss_fft_cpx* image, kiss_fft_cpx* kspace, kiss_fft_cpx* filt_image, kiss_fft_cpx* filt_kspace, kiss_fft_cpx* cs_image, kiss_fft_cpx* cs_kspace, dataset* ds) {
    int x,y,z,ipos,ds_pos,nx,ny,nz;
    double xs,ys,zs;
    double s,d;
    xs = (double)ds->k_xdim/(double)p->xdim/2.0;
    ys = (double)ds->k_ydim/(double)p->ydim/2.0;
    zs = (double)ds->k_zdim/(double)p->zdim/2.0;
    /*fprintf(stderr, "%i %i ", p->sequence, p->n_params);
    for(int i=0;i<p->n_params;i++) {
        fprintf(stderr, "%f ", p->s_params[i]);
    }
    fprintf(stderr, "%i %i %i\n", p->xdim, p->ydim, p->zdim);*/
    for(z = 0; z<p->zdim; z++) {
        for(y = 0;y<p->ydim; y++) {
            for(x = 0;x<p->xdim; x++) {
                ipos = z*p->xdim*p->ydim + y*p->xdim + x;
                d = 0; s = 0;
                nz = z*ds->k_zdim/p->zdim;
                ny = y*ds->k_ydim/p->ydim;
                nx = x*ds->k_xdim/p->xdim;
                switch (p->nearest) {
                    case Nearest:
                        ds_pos = round(nz)*ds->k_xdim*ds->k_ydim+round(ny)*ds->k_xdim+round(nx);
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
                        d = 0; s=0;
                        for(int xi=floor(nx-xs);xi<ceil(nx+xs);xi+=1) {
                            for(int yi=floor(ny-ys);yi<ceil(ny+ys);yi+=1) {
                                for(int zi=floor(nz-zs);zi<ceil(nz+zs);zi+=1) {
                                    ds_pos = zi*ds->k_xdim*ds->k_ydim+yi*ds->k_xdim+xi;
                                    if(ds_pos >= 0 && ds_pos < ds->k_xdim*ds->k_ydim*ds->k_zdim && xi>=0 && xi<ds->k_xdim && yi>=0 && yi<ds->k_ydim && zi>=0 && zi<ds->k_zdim) {
                                        s += simVoxel(p, ds_pos, ds);
                                        d += 1.0;
                                    }
                                }
                            }
                        }
                        s /= d;
                        break;
                }
                ASSIGN(image[ipos], (kiss_fft_scalar)s, 0);
            }
        }
    }
    //free_dataset(ds);

    if(p->use_fft3) {
        fft3d(image, kspace, p->xdim, p->ydim, p->zdim);
    } else {
        for(z=0;z<p->zdim;z++) {
            fft2d(&image[z*p->xdim*p->ydim], &kspace[z*p->xdim*p->ydim], p->xdim, p->ydim);
        }
    }

    if(p->use_cs && cs_image != NULL) {
        apply_cs_filter(kspace, filt_kspace, p->cs_params);
        if(p->use_fft3) {
            ifft3d(filt_kspace, filt_image, p->xdim, p->ydim, p->zdim);
        } else {
            for(z=0;z<p->zdim;z++) {
                ifft2d(&filt_kspace[z*p->xdim*p->ydim], &filt_image[z*p->xdim*p->ydim], p->xdim, p->ydim);
            }
        }
        for(int i=0;i<p->xdim*p->ydim*p->zdim;i++) {
            cs_kspace[i] = filt_kspace[i];
            ASSIGN(cs_image[i], 0, 0);
        }
        p->cs_params->zdim = 1;
        for(z=0;z<p->zdim;z++) {
            if (p->cs_params->callback != 0) {
                ((cs_callback*)p->cs_params->callback)(z);
            } else {
                fprintf(stderr, "%d\n", z);
            }
            compressed_sensing(&cs_kspace[z*p->xdim*p->ydim], &cs_image[z*p->xdim*p->ydim], p->cs_params);
        }
    }
}