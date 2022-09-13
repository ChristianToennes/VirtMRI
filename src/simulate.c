#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include "fft.h"
#include "fft2.h"
#include "fft3.h"
#include "tinycs.h"
#include "noise.h"
#include "stdbool.h"

#include "_kiss_fft_guts.h"
#include "params.h"

#define NA_SCALE 140

void free_dataset(struct Dataset *ds) {
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

static inline double simVoxel(struct Params *p, int pos, struct Dataset *ds) {
    double te, tr, ti, fa, tau1, tau2, te_end, te_step, e_tr_t1, e_tr_t2, cfa, sfa, s, m, b;
    double pd, t1, t2, t2s, t2f, ex_frac;
    s = 0;
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
            s = fabs(pd*sfa*(1.0-e_tr_t1)/(1.0-(e_tr_t1+e_tr_t2)*cfa-e_tr_t1*e_tr_t2)*exp(-te/t2) );
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
            s = fabs((na_vol-ex_frac)*pd* (1-exp(-tr/t1)) * (0.6*exp(-te/t2f) + 0.4*exp(-te/t2s)) + ex_frac*pd* (1-exp(-tr/t1))*exp(-te/na_t2fr))/NA_SCALE;
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
            for(;te<=te_end;te+=te_step) {
                s += fabs(pd*(exp(-(te+tau1)/t2s) * exp(-(te+tau1)/t2f))*sfa) / NA_SCALE;
            }
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
            for(;te<=te_end;te+=te_step) {
                s += fabs(pd*( (exp(-te/t2s) - exp(-te/t2f)) * (exp(-tau1/t2s)-exp(-tau1/t2f)) * exp(-tau2/t2s) )) / NA_SCALE;
            }
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
            for(;te<=te_end;te+=te_step) {
                e_tr_t1 += fabs(pd*( (exp(-te/t2s)-exp(-te/t2f)) * (exp(-tau1/t2s)-exp(-tau1/t2f)) * exp(-tau2/t2s) )) / NA_SCALE;
            }
            fa = p->s_params[6] * M_PI / 180;
            tau1 = p->s_params[7];
            te = p->s_params[8];
            te_end = p->s_params[9];
            te_step = p->s_params[10];
            sfa = sin(fa);
            e_tr_t2 = 0;
            for(;te<=te_end;te+=te_step) {
                e_tr_t2 += fabs(pd*(exp(-(te+tau1)/t2s) * exp(-(te+tau1)/t2f))*sfa) / NA_SCALE;
            }
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
            s = fabs(pd*(exp(-te/t2s)-exp(-te/t2f)) * (exp(-tau1/t2s)-exp(-tau1/t2f)) * exp(-tau2/t2s) * (sfa*sfa*sfa*sfa*sfa));
            break;
    }
    return s;
}

void simulate(struct Params *p, kiss_fft_cpx *image, kiss_fft_cpx *kspace, kiss_fft_cpx *filt_image, kiss_fft_cpx *filt_kspace, kiss_fft_cpx *cs_image, kiss_fft_cpx *cs_kspace, struct Dataset *ds) {
    int x,y,z,ipos,ds_pos,nx,ny,nz;
    double xs,ys,zs;
    double s,d;
    xs = (double)ds->k_xdim/(double)p->xdim/2.0;
    ys = (double)ds->k_ydim/(double)p->ydim/2.0;
    zs = (double)ds->k_zdim/(double)p->zdim/2.0;
    
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

    addImageNoise(image, p);

    if(p->use_fft3) {
        fft3d(image, kspace, p->ydim, p->zdim, p->xdim);
    } else {
        for(z=0;z<p->zdim;z++) {
            fft2d(&image[z*p->xdim*p->ydim], &kspace[z*p->xdim*p->ydim], p->xdim, p->ydim);
        }
    }

    if(addKSpaceNoise(kspace, p)) {
        if(p->use_fft3) {
            ifft3d(kspace, image, p->ydim, p->zdim, p->xdim);
        } else {
            for(z=0;z<p->zdim;z++) {
                ifft2d(&kspace[z*p->xdim*p->ydim], &image[z*p->xdim*p->ydim], p->xdim, p->ydim);
            }
        }
    }

    if(p->use_cs && cs_image != NULL) {
        apply_cs_filter(kspace, filt_kspace, p->cs_params);
        if(p->use_fft3) {
            ifft3d(filt_kspace, filt_image, p->ydim, p->zdim, p->xdim);
        } else {
            for(z=0;z<p->zdim;z++) {
                ifft2d(&filt_kspace[z*p->xdim*p->ydim], &filt_image[z*p->xdim*p->ydim], p->xdim, p->ydim);
            }
        }
        for(int i=0;i<p->xdim*p->ydim*p->zdim;i++) {
            cs_kspace[i] = filt_kspace[i];
            ASSIGN(cs_image[i], 0, 0);
        }
        if(p->use_fft3) {
            compressed_sensing(cs_kspace, cs_image, p->cs_params);
        } else {
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
}