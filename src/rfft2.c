#include <stdint.h>
#include <stdio.h>

#include "kiss_fftndr.h"

void rfft2d(float *in, float *out, int m, int n) {
   int dims[] = {n, m};
   kiss_fftndr_cfg cfg;

   cfg = kiss_fftndr_alloc(dims, 2, 0, NULL, NULL);
   kiss_fftndr(cfg, in, (kiss_fft_cpx *) out);
   free(cfg);
}

void irfft2d(float *in, float *out, int m, int n) {
   int dims[] = {n, m};
   kiss_fftndr_cfg cfg;

   cfg = kiss_fftndr_alloc(dims, 2, 1, NULL, NULL);
   kiss_fftndri(cfg, (const kiss_fft_cpx *) in, out);
   free(cfg);
}
