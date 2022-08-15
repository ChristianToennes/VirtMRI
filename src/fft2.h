#ifndef __fft2_h
#define __fft2_h

#include "kiss_fftnd.h"

void fft2d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n);

void ifft2d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n);

#endif