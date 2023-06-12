#ifndef __fft2_h
#define __fft2_h

#include "kiss_fftnd.h"

void kfft2d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n);

void kifft2d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n);

#endif