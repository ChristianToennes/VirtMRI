#include "tinycs.h"

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include "stdbool.h"
#include "kissfft/_kiss_fft_guts.h"

void print_stats(char* name, kiss_fft_cpx* data, int length) {
    fprintf(stderr, "%s", name);
    int nans = 0;
    double sum = 0;
    kiss_fft_scalar min, max;
    min = REAL(data[0]);
    max = REAL(data[0]);
    for(int i=0;i<length;i++) {
        if(REAL(data[i])!=REAL(data[i]) || IMAG(data[i])!=IMAG(data[i])) {
            nans++;
        }
        if(REAL(data[i])>max) { max = REAL(data[i]); }
        if(REAL(data[i])<min) { min = REAL(data[i]); }
        sum += REAL(data[i])+IMAG(data[i]);
    }
    fprintf(stderr, " %i %f %f %f\n", nans, min, max, sum);
}

void compressed_sensing(kiss_fft_cpx *f_data, kiss_fft_cpx *out, struct Params *params) {
    int xdim = params->xdim;
    int ydim = params->ydim;
    int zdim = params->zdim;
    if(!params->use_fft3) {
        zdim = 1;
    }
    int k_xdim = xdim;
    int k_ydim = ydim;
    int k_zdim = zdim;
    int ninner = params->cs_params->ninner;
    int nbreg = params->cs_params->nbreg;
    double lambda = params->cs_params->lambda;
    double lambda2 = params->cs_params->lambda2;
    double mu = params->cs_params->mu;
    double gam = params->cs_params->gam;
    int data_ndims = zdim == 1 ? 2 : 3;
    int data_length = xdim*ydim*zdim;
    char* f_mask = (char*)malloc(data_length);
    for(int i=0;i<data_length;i++) {
        f_mask[i] = (REAL(f_data[i])+IMAG(f_data[i])) != 0;
    }
    //% normalize the data so that standard parameter values work
    double norm_factor;
    norm_factor = 0;
    kiss_fft_cpx tmp;
    for(int i=0;i<data_length;i++) {
        norm_factor += REAL(f_data[i])*REAL(f_data[i])+IMAG(f_data[i])*IMAG(f_data[i]);
    }
    //fprintf(stderr, "sum %f\n", norm_factor);
    norm_factor = 1.0 / (sqrt(norm_factor) / (double)ydim);
    //fprintf(stderr, "sqrt %f\n", norm_factor);
    for(int i=0;i<data_length;i++) {
        C_MULBYSCALAR(f_data[i], norm_factor);
    }
    //print_stats("f_data", f_data, data_length);
    //% Reserve memory for the auxillary variables
    kiss_fft_cpx* f_data0 = (kiss_fft_cpx*)malloc(data_length*sizeof(kiss_fft_cpx));
    for(int i=0;i<data_length;i++) {f_data0[i] = f_data[i];}
    kiss_fft_cpx* img = (kiss_fft_cpx*)malloc(data_length*sizeof(kiss_fft_cpx));
    kiss_fft_cpx* k_img = (kiss_fft_cpx*)malloc(data_length*sizeof(kiss_fft_cpx));
    kiss_fft_cpx* X = (kiss_fft_cpx*)malloc(data_length*sizeof(kiss_fft_cpx)*3);
    kiss_fft_cpx* B = (kiss_fft_cpx*)malloc(data_length*sizeof(kiss_fft_cpx)*3);
    for(int i=0;i<data_length*3;i++) {
        ASSIGN(X[i], 0, 0);
        ASSIGN(B[i], 0, 0);
    }
    //fprintf(stderr, "assign\n");
    kiss_fft_cpx* i_murf = (kiss_fft_cpx*)malloc(data_length*sizeof(kiss_fft_cpx));
    kiss_fft_cpx* i_rhs = (kiss_fft_cpx*)malloc(data_length*sizeof(kiss_fft_cpx));
    kiss_fft_cpx* f_rhs = (kiss_fft_cpx*)malloc(data_length*sizeof(kiss_fft_cpx));
    //% Build Kernels
    double scale = sqrt(data_length);
    scale = 1;
    kiss_fft_cpx* f_uker = (kiss_fft_cpx*)malloc(data_length*sizeof(kiss_fft_cpx));
    //fprintf(stderr, "%i %i %i %i\n", f_uker, i_rhs, img, k_img);
    for(int i=0;i<data_length;i++) {
        ASSIGN(img[i], 0, 0);
        ASSIGN(k_img[i], 0, 0);
        ASSIGN(i_murf[i], 0, 0);
        ASSIGN(i_rhs[i], 0, 0);
        ASSIGN(f_rhs[i], 0, 0);
        ASSIGN(f_uker[i], 0, 0);
    }
    if (data_ndims == 2) {
      ASSIGN(f_uker[0], 4, 0);
      ASSIGN(f_uker[1], -1, 0);
      ASSIGN(f_uker[xdim], -1, 0);
      ASSIGN(f_uker[(xdim-1)], -1, 0);
      ASSIGN(f_uker[xdim*(ydim-1)], -1, 0);
    } else {// data_ndims == 3
      ASSIGN(f_uker[0], 8, 0);
      ASSIGN(f_uker[1], -1, 0);
      ASSIGN(f_uker[xdim], -1, 0);
      ASSIGN(f_uker[xdim*ydim], -1, 0);
      ASSIGN(f_uker[(xdim-1)], -1, 0);
      ASSIGN(f_uker[xdim*(ydim-1)], -1, 0);
      ASSIGN(f_uker[xdim*ydim*(zdim-1)], -1, 0);
    }
    //print_stats("uker in", f_uker, data_length);
    data_ndims==2 ? fft2d(f_uker, f_uker, xdim, ydim) : fft3d(f_uker, f_uker, ydim, zdim, xdim);
    //print_stats("uker mid", f_uker);
    //console.log(mu, lambda, gam);
    //fprintf(stderr, "fft\n");
    kiss_fft_scalar real, imag, divis;
    for(int i=0;i<data_length;i++) {
        real = (mu*f_mask[i]+lambda*REAL(f_uker[i])+gam);
        imag = lambda*IMAG(f_uker[i]);
        divis = real*real+imag*imag;
        ASSIGN(f_uker[i], real/divis, -imag/divis);
    }
    //print_stats("uker out", f_uker, data_length);
    //%  Do the reconstruction
    kiss_fft_cpx d;
    int i, x, y, z;
    for(int outer = 0;outer < nbreg; outer++) {
      if (!(data_ndims==2) && params->cs_params->callback != 0) {
        ((cs_callback*)params->cs_params->callback)(outer+1);
      }
      for(int i=0;i<data_length;i++) {
        ASSIGN(i_rhs[i], REAL(f_data[i]), IMAG(f_data[i]));
        C_MULBYSCALAR(i_rhs[i], mu*f_mask[i]);
      }
      //print_stats("murf in", i_rhs, data_length);
      data_ndims==2? ifft2d(i_rhs, i_murf, k_xdim, k_ydim) : ifft3d(i_rhs, i_murf, k_ydim, k_zdim, k_xdim);
      //print_stats("murf out", i_murf, data_length);
      for(int inner = 0;inner < ninner;inner++) {
        //% update u
        for(x=0;x<xdim;x++) {
            for(y=0;y<ydim;y++) {
                for(z=0;z<zdim;z++) {
                    i = 3*(x + y*xdim + z*xdim*ydim);
                    ASSIGN(d, 0, 0);
                    if(x<xdim-1) {
                        C_ADDTO(d, X[i]);
                        C_SUBFROM(d, B[i]);
                        C_SUBFROM(d, X[3*(x+1 + y*xdim + z*xdim*ydim)]);
                        C_ADDTO(d, B[3*(x+1 + y*xdim + z*xdim*ydim)]);
                    } else {
                        C_ADDTO(d, X[i]);
                        C_SUBFROM(d, B[i]);
                        C_SUBFROM(d, X[3*(x+1 + y*xdim + z*xdim*ydim)]);
                        C_ADDTO(d, B[3*(x+1 + y*xdim + z*xdim*ydim)]);
                    }
                    if(y<ydim-1) {
                        C_ADDTO(d, X[i+1]);
                        C_SUBFROM(d, B[i+1]);
                        C_SUBFROM(d, X[3*(x + (y+1)*xdim + z*xdim*ydim)+1]);
                        C_ADDTO(d, B[3*(x + (y+1)*xdim + z*xdim*ydim)+1]);
                    } else {
                        C_ADDTO(d, X[i+1]);
                        C_SUBFROM(d, B[i+1]);
                        C_SUBFROM(d, X[3*(x + z*xdim*ydim)+1]);
                        C_ADDTO(d, B[3*(x + z*xdim*ydim)+1]);
                    }
                    if(z<zdim-1) {
                        C_ADDTO(d, X[i+2]);
                        C_SUBFROM(d, B[i+2]);
                        C_SUBFROM(d, X[3*(x + y*xdim + (z+1)*xdim*ydim)+2]);
                        C_ADDTO(d, B[3*(x + y*xdim + (z+1)*xdim*ydim)+2]);
                    } else {
                        C_ADDTO(d, X[i+2]);
                        C_SUBFROM(d, B[i+2]);
                        C_SUBFROM(d, X[3*(x + y*xdim)+2]);
                        C_ADDTO(d, B[3*(x + y*xdim)+2]);
                    }
                    i = x + y*xdim + z*xdim*ydim;
                    //fprintf(stderr, "d %f+%fi", REAL(d), IMAG(d));
                    ASSIGN(i_rhs[i], REAL(i_murf[i]), IMAG(i_murf[i]));
                    //C_MULBYSCALAR(i_rhs[i], scale);
                    C_MULBYSCALAR(d, lambda);
                    C_ADDTO(i_rhs[i], d);
                    ASSIGN(tmp, REAL(img[i]), IMAG(img[i]));
                    C_MULBYSCALAR(tmp, gam);
                    C_ADDTO(i_rhs[i], tmp);
                }
            }
        }
        //fprintf(stderr, "%f\n", scale);
        //print_stats("rhs in", i_rhs, data_length);
        data_ndims==2? fft2d(i_rhs, f_rhs, k_xdim, k_ydim) : fft3d(i_rhs, f_rhs, k_ydim, k_zdim, k_xdim);
        //print_stats("rhs out", f_rhs, data_length);
        for(i=0;i<data_length;i++) {
            C_MUL(tmp, f_rhs[i], f_uker[i]);
            ASSIGN(f_rhs[i], REAL(tmp), IMAG(tmp));
        }
        //print_stats("img_in", f_rhs, data_length);
        data_ndims==2? ifft2d(f_rhs, img, k_xdim, k_ydim) : ifft3d(f_rhs, img, k_ydim, k_zdim, k_xdim);
        //print_stats("img", img, data_length);
        //% update x and y
        for(x=0;x<xdim;x++) {
            for(y=0;y<ydim;y++) {
                for(z=0;z<zdim;z++) {
                    i = 3*(x + y*xdim + z*xdim*ydim);
                    if(x>0) {
                        C_ADDTO(B[i], img[x + y*xdim + z*xdim*ydim]);
                        C_SUBFROM(B[i], img[x-1 + y*xdim + z*xdim*ydim]);
                    } else {
                        C_ADDTO(B[i], img[x + y*xdim + z*xdim*ydim]);
                        C_SUBFROM(B[i], img[xdim-1 + y*xdim + z*xdim*ydim]);
                    }
                    if(y>0) {
                        C_ADDTO(B[i+1], img[x + y*xdim + z*xdim*ydim]);
                        C_SUBFROM(B[i+1], img[x + (y-1)*xdim + z*xdim*ydim]);
                    } else {
                        C_ADDTO(B[i+1], img[x + y*xdim + z*xdim*ydim]);
                        C_SUBFROM(B[i+1], img[x + (ydim-1)*xdim + z*xdim*ydim]);
                    }
                    if(z>0) {
                        C_ADDTO(B[i+2], img[x + y*xdim + z*xdim*ydim]);
                        C_SUBFROM(B[i+2], img[x + y*xdim + (z-1)*xdim*ydim]);
                    } else {
                        C_ADDTO(B[i+2], img[x + y*xdim + z*xdim*ydim]);
                        C_SUBFROM(B[i+2], img[x + y*xdim + (zdim-1)*xdim*ydim]);
                    }
                }
            }
        }
        for(i=0;i<data_length*3;i+=3) {
            double s = sqrt(REAL(B[i])*REAL(B[i])+IMAG(B[i])*IMAG(B[i])+REAL(B[i+1])*REAL(B[i+1])+IMAG(B[i+1])*IMAG(B[i+1])+REAL(B[i+2])*REAL(B[i+2])+IMAG(B[i+2])*IMAG(B[i+2]));
            double ss = s - lambda2;
            ss = ss > 0 ? ss : 0.0;
            s = s + (s < lambda2 ? 1.0 : 0.0);
            ss = ss / s;
            ASSIGN(X[i], REAL(B[i]), IMAG(B[i]));
            C_MULBYSCALAR(X[i], ss);
            ASSIGN(X[i+1], REAL(B[i+1]), IMAG(B[i+1]));
            C_MULBYSCALAR(X[i+1], ss);
            ASSIGN(X[i+2], REAL(B[i+2]), IMAG(B[i+2]));
            C_MULBYSCALAR(X[i+2], ss);
            //% update bregman parameters
            C_SUBFROM(B[i], X[i]);
            C_SUBFROM(B[i+1], X[i+1]);
            C_SUBFROM(B[i+2], X[i+2]);
        }
      }
      data_ndims==2? fft2d(img, k_img, k_xdim, k_ydim) : fft3d(img, k_img, k_ydim, k_zdim, k_xdim);
      //print_stats("k_img", k_img, data_length);
      for(i=0;i<data_length;i++) {
        C_ADDTO(f_data[i], f_data0[i]);
        ASSIGN(tmp, REAL(k_img[i]), IMAG(k_img[i]));
        C_MULBYSCALAR(tmp, f_mask[i] / scale);
        C_SUBFROM(f_data[i], tmp); 
      }
    }
    //% undo the normalization so that results are scaled properly
    //print_stats("img in", img, data_length);
    //fprintf(stderr, "make img %f %f\n", norm_factor, scale);
    for(i=0;i<data_length;i++) {
        ASSIGN(out[i], sqrt(REAL(img[i])*REAL(img[i])+IMAG(img[i])*IMAG(img[i])) / norm_factor / scale, 0);
        ASSIGN(f_data[i], REAL(k_img[i]), IMAG(k_img[i]));
        //out[i].i = 0;
        //out[i] = img[i];
    }
    //fprintf(stderr, "%i\n", data_length);
    free(f_data0);
    free(img);
    free(k_img);
    free(X);
    free(B);
    free(i_murf);
    free(i_rhs);
    free(f_rhs);
    free(f_uker);
    //free(f_mask);
    //print_stats("img", out, data_length);
}
