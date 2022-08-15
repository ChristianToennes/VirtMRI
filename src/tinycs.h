#ifndef __tinycs_h
#define __tinycs_h

#include "kissfft/kiss_fft.h"

void fft2d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n);
void ifft2d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n);
void fft3d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n, int l);
void ifft3d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n, int l);

typedef struct {
    int xdim;
    int ydim;
    int zdim;
    int ninner;
    int nbreg;
    double lambda;
    double lambda2;
    double mu;
    double gam;
} cs_params;

cs_params* make_cs_params(int xdim, int ydim, int zdim, int ninner, int nbreg, double lambda, double lambda2, double mu, double gam);

void print_stats(char* name, kiss_fft_cpx* data, int length);

void compressed_sensing(kiss_fft_cpx *f_data, kiss_fft_cpx *out, cs_params *params);

#endif