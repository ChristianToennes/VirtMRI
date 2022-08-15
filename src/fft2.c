#include "fft2.h"

#include <stdint.h>
#include <stdio.h>

#include "_kiss_fft_guts.h"

void fft2d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n) {
   int dims[] = {n, m};
   kiss_fftnd_cfg cfg;

   cfg = kiss_fftnd_alloc(dims, 2, 0, NULL, NULL);
   kiss_fftnd(cfg, in, out);
   free(cfg);
}

void ifft2d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n) {
   int dims[] = {n, m};
   kiss_fftnd_cfg cfg;

   cfg = kiss_fftnd_alloc(dims, 2, 1, NULL, NULL);
   kiss_fftnd(cfg, in, out);
   for(int i=0;i<m*n;i++) {
      C_MULBYSCALAR(out[i], 1.0/(m*n));
   }
   free(cfg);
}
