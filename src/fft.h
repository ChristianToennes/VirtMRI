#ifndef __fft_h
#define __fft_h
#include "kiss_fft.h"

void kfft(kiss_fft_cpx *in, kiss_fft_cpx *out, int m);

void kifft(kiss_fft_cpx *in, kiss_fft_cpx *out, int m);

void kfft_stride(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int nl);

#endif