#ifndef __fft3_h
#define __fft3_h

#include "kiss_fftnd.h"

void __fft3d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n, int l);

void __ifft3d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n, int l);

void fft3d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n, int l);

void ifft3d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n, int l);

#endif