#ifndef __fft_h
#define __fft_h
#include "kiss_fft.h"

void fft(kiss_fft_cpx *in, kiss_fft_cpx *out, int m);

void ifft(kiss_fft_cpx *in, kiss_fft_cpx *out, int m);

void fft_stride(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int nl);

#endif