#include "tinycs.h"

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include "stdbool.h"

#include "num/fft.h"

void print_stats(char* name, complex float* data, int length) {
    fprintf(stderr, "%s", name);
    int nans = 0;
    double sum = 0;
    float min, max;
    min = creal(data[0]);
    max = creal(data[0]);
    for(int i=0;i<length;i++) {
        if(creal(data[i])!=creal(data[i]) || cimag(data[i])!=cimag(data[i])) {
            nans++;
        }
        if(creal(data[i])>max) { max = creal(data[i]); }
        if(creal(data[i])<min) { min = creal(data[i]); }
        sum += creal(data[i])+cimag(data[i]);
    }
    fprintf(stderr, " %i %f %f %f\n", nans, min, max, sum);
}

void compressed_sensing(complex float *f_data, complex float *out, struct Params *params) {
    int xdim = params->ixdim;
    int ydim = params->iydim;
    int zdim = params->izdim;
    if(!params->use_fft3) {
        zdim = 1;
    }
    int k_xdim = xdim;
    int k_ydim = ydim;
    int k_zdim = zdim;
    long dims3[] = {params->izdim, params->iydim, params->ixdim};
    long dims2[] = {params->iydim, params->ixdim};
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
        f_mask[i] = (creal(f_data[i])+cimag(f_data[i])) != 0;
    }
    //% normalize the data so that standard parameter values work
    double norm_factor;
    norm_factor = 0;
    complex float tmp;
    for(int i=0;i<data_length;i++) {
        norm_factor += creal(f_data[i])*creal(f_data[i])+cimag(f_data[i])*cimag(f_data[i]);
    }
    //fprintf(stderr, "sum %f\n", norm_factor);
    norm_factor = 1.0 / (sqrt(norm_factor) / (double)ydim);
    //fprintf(stderr, "sqrt %f\n", norm_factor);
    for(int i=0;i<data_length;i++) {
        f_data[i] *= norm_factor;
    }
    //print_stats("f_data", f_data, data_length);
    //% Reserve memory for the auxillary variables
    complex float* f_data0 = (complex float*)malloc(data_length*sizeof(complex float));
    for(int i=0;i<data_length;i++) {f_data0[i] = f_data[i];}
    complex float* img = (complex float*)malloc(data_length*sizeof(complex float));
    complex float* k_img = (complex float*)malloc(data_length*sizeof(complex float));
    complex float* X = (complex float*)malloc(data_length*sizeof(complex float)*3);
    complex float* B = (complex float*)malloc(data_length*sizeof(complex float)*3);
    for(int i=0;i<data_length*3;i++) {
        X[i] = 0+0*I;
        B[i] = 0+0*I;
    }
    //fprintf(stderr, "assign\n");
    complex float* i_murf = (complex float*)malloc(data_length*sizeof(complex float));
    complex float* i_rhs = (complex float*)malloc(data_length*sizeof(complex float));
    complex float* f_rhs = (complex float*)malloc(data_length*sizeof(complex float));
    //% Build Kernels
    double scale = sqrt(data_length);
    scale = 1;
    complex float* f_uker = (complex float*)malloc(data_length*sizeof(complex float));
    //fprintf(stderr, "%i %i %i %i\n", f_uker, i_rhs, img, k_img);
    for(int i=0;i<data_length;i++) {
        img[i] = 0+0*I;
        k_img[i] = 0+0*I;
        i_murf[i] = 0+0*I;
        i_rhs[i] = 0+0*I;
        f_rhs[i] = 0+0*I;
        f_uker[i] = 0+0*I;
    }
    if (data_ndims == 2) {
      f_uker[0] = 4+0*I;
      f_uker[1] = -1+0*I;
      f_uker[xdim] = -1+0*I;
      f_uker[(xdim-1)] = -1+0*I;
      f_uker[xdim*(ydim-1)] = -1+0*I;
    } else {// data_ndims == 3
      f_uker[0] = 8+0*I;
      f_uker[1] = -1+0*I;
      f_uker[xdim] = -1+0*I;
      f_uker[xdim*ydim] = -1+0*I;
      f_uker[(xdim-1)] = -1+0*I;
      f_uker[xdim*(ydim-1)] = -1+0*I;
      f_uker[xdim*ydim*(zdim-1)] = -1+0*I;
    }
    //print_stats("uker in", f_uker, data_length);
    data_ndims==2 ? fft(2, dims2, 3, f_uker, f_uker) : fft(3, dims3, 7, f_uker, f_uker);
    //data_ndims==2 ? kfft2d(f_uker, f_uker, xdim, ydim) : kfft3d(f_uker, f_uker, ydim, zdim, xdim);
    //data_ndims==2 ? kfft2d(f_uker, f_uker, xdim, ydim) : kfft3d(f_uker, f_uker, ydim, zdim, xdim);
    //print_stats("uker mid", f_uker);
    //console.log(mu, lambda, gam);
    //fprintf(stderr, "fft\n");
    float creal, cimag, divis;
    for(int i=0;i<data_length;i++) {
        creal = (mu*f_mask[i]+lambda*creal(f_uker[i])+gam);
        cimag = lambda*cimag(f_uker[i]);
        divis = creal*creal+cimag*cimag;
        f_uker[i] = creal/divis+-cimag/divis*I;
    }
    //print_stats("uker out", f_uker, data_length);
    //%  Do the reconstruction
    complex float d;
    int i, x, y, z;
    for(int outer = 0;outer < nbreg; outer++) {
      if (!(data_ndims==2) && params->cs_params->callback != 0) {
        ((cs_callback*)params->cs_params->callback)(outer+1);
      }
      for(int i=0;i<data_length;i++) {
        i_rhs[i] = creal(f_data[i])+cimag(f_data[i])*I;
        i_rhs[i] *= mu*f_mask[i];
      }
      //print_stats("murf in", i_rhs, data_length);
      data_ndims==2 ? fft(2, dims2, 3, i_murf, i_rhs) : fft(3, dims3, 7, i_murf, i_rhs);
      //data_ndims==2? kifft2d(i_rhs, i_murf, k_xdim, k_ydim) : kifft3d(i_rhs, i_murf, k_ydim, k_zdim, k_xdim);
      //print_stats("murf out", i_murf, data_length);
      for(int inner = 0;inner < ninner;inner++) {
        //% update u
        for(x=0;x<xdim;x++) {
            for(y=0;y<ydim;y++) {
                for(z=0;z<zdim;z++) {
                    i = 3*(x + y*xdim + z*xdim*ydim);
                    d = 0+0*I;
                    if(x<xdim-1) {
                        d += X[i] ;
                        d -= B[i] ;
                        d -= X[3*(x+1 + y*xdim + z*xdim*ydim)] ;
                        d += B[3*(x+1 + y*xdim + z*xdim*ydim)] ;
                    } else {
                        d += X[i] ;
                        d -= B[i] ;
                        d -= X[3*(x+1 + y*xdim + z*xdim*ydim)] ;
                        d += B[3*(x+1 + y*xdim + z*xdim*ydim)] ;
                    }
                    if(y<ydim-1) {
                        d += X[i+1] ;
                        d -= B[i+1] ;
                        d -= X[3*(x + (y+1)*xdim + z*xdim*ydim)+1] ;
                        d += B[3*(x + (y+1)*xdim + z*xdim*ydim)+1] ;
                    } else {
                        d += X[i+1] ;
                        d -= B[i+1] ;
                        d -= X[3*(x + z*xdim*ydim)+1] ;
                        d += B[3*(x + z*xdim*ydim)+1] ;
                    }
                    if(z<zdim-1) {
                        d += X[i+2] ;
                        d -= B[i+2] ;
                        d -= X[3*(x + y*xdim + (z+1)*xdim*ydim)+2] ;
                        d += B[3*(x + y*xdim + (z+1)*xdim*ydim)+2] ;
                    } else {
                        d += X[i+2] ;
                        d -= B[i+2] ;
                        d -= X[3*(x + y*xdim)+2] ;
                        d += B[3*(x + y*xdim)+2] ;
                    }
                    i = x + y*xdim + z*xdim*ydim;
                    //fprintf(stderr, "d %f+%fi", creal(d), cimag(d));
                    i_rhs[i] = creal(i_murf[i])+cimag(i_murf[i])*I;
                    //C_MULBYSCALAR(i_rhs[i], scale);
                    d *= lambda;
                    i_rhs[i] += d ;
                    tmp = creal(img[i])+cimag(img[i])*I;
                    tmp *= gam;
                    i_rhs[i] += tmp ;
                }
            }
        }
        //fprintf(stderr, "%f\n", scale);
        //print_stats("rhs in", i_rhs, data_length);
        data_ndims==2 ? fft(2, dims2, 3, f_rhs, i_rhs) : fft(3, dims3, 7, f_rhs, i_rhs);
        //data_ndims==2? kfft2d(i_rhs, f_rhs, k_xdim, k_ydim) : kfft3d(i_rhs, f_rhs, k_ydim, k_zdim, k_xdim);
        //print_stats("rhs out", f_rhs, data_length);
        for(i=0;i<data_length;i++) {
            tmp = f_rhs[i] * f_uker[i];
            f_rhs[i] = creal(tmp)+cimag(tmp)*I;
        }
        //print_stats("img_in", f_rhs, data_length);
        data_ndims==2 ? fft(2, dims2, 3, img, f_rhs) : fft(3, dims3, 7, img, f_rhs);
        //data_ndims==2? kifft2d(f_rhs, img, k_xdim, k_ydim) : kifft3d(f_rhs, img, k_ydim, k_zdim, k_xdim);
        //print_stats("img", img, data_length);
        //% update x and y
        for(x=0;x<xdim;x++) {
            for(y=0;y<ydim;y++) {
                for(z=0;z<zdim;z++) {
                    i = 3*(x + y*xdim + z*xdim*ydim);
                    if(x>0) {
                        B[i] += img[x + y*xdim + z*xdim*ydim] ;
                        B[i] -= img[x-1 + y*xdim + z*xdim*ydim] ;
                    } else {
                        B[i] += img[x + y*xdim + z*xdim*ydim] ;
                        B[i] -= img[xdim-1 + y*xdim + z*xdim*ydim] ;
                    }
                    if(y>0) {
                        B[i+1] += img[x + y*xdim + z*xdim*ydim] ;
                        B[i+1] -= img[x + (y-1)*xdim + z*xdim*ydim] ;
                    } else {
                        B[i+1] += img[x + y*xdim + z*xdim*ydim] ;
                        B[i+1] -= img[x + (ydim-1)*xdim + z*xdim*ydim] ;
                    }
                    if(z>0) {
                        B[i+2] += img[x + y*xdim + z*xdim*ydim] ;
                        B[i+2] -= img[x + y*xdim + (z-1)*xdim*ydim] ;
                    } else {
                        B[i+2] += img[x + y*xdim + z*xdim*ydim] ;
                        B[i+2] -= img[x + y*xdim + (zdim-1)*xdim*ydim] ;
                    }
                }
            }
        }
        for(i=0;i<data_length*3;i+=3) {
            double s = sqrt(creal(B[i])*creal(B[i])+cimag(B[i])*cimag(B[i])+creal(B[i+1])*creal(B[i+1])+cimag(B[i+1])*cimag(B[i+1])+creal(B[i+2])*creal(B[i+2])+cimag(B[i+2])*cimag(B[i+2]));
            double ss = s - lambda2;
            ss = ss > 0 ? ss : 0.0;
            s = s + (s < lambda2 ? 1.0 : 0.0);
            ss = ss / s;
            X[i] = creal(B[i])+cimag(B[i])*I;
            X[i] *= ss;
            X[i+1] = creal(B[i+1])+cimag(B[i+1])*I;
            X[i+1] *= ss;
            X[i+2] = creal(B[i+2])+cimag(B[i+2])*I;
            X[i+2] *= ss;
            //% update bregman parameters
            B[i] -= X[i] ;
            B[i+1] -= X[i+1] ;
            B[i+2] -= X[i+2] ;
        }
      }
      data_ndims==2 ? fft(2, dims2, 3, k_img, img) : fft(3, dims3, 7, k_img, img);
      //data_ndims==2? kfft2d(img, k_img, k_xdim, k_ydim) : kfft3d(img, k_img, k_ydim, k_zdim, k_xdim);
      //print_stats("k_img", k_img, data_length);
      for(i=0;i<data_length;i++) {
        f_data[i] += f_data0[i] ;
        tmp = creal(k_img[i])+cimag(k_img[i])*I;
        tmp *= f_mask[i] / scale;
        f_data[i] -= tmp ; 
      }
    }
    //% undo the normalization so that results are scaled properly
    //print_stats("img in", img, data_length);
    //fprintf(stderr, "make img %f %f\n", norm_factor, scale);
    for(i=0;i<data_length;i++) {
        out[i] = sqrt(creal(img[i])*creal(img[i])+cimag(img[i])*cimag(img[i])) / norm_factor / scale+0*I;
        f_data[i] = creal(k_img[i])+cimag(k_img[i])*I;
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
