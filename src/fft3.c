#include "fft3.h"

#include <stdint.h>
#include <stdio.h>

#include "_kiss_fft_guts.h"

void fft3d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n, int l) {
   int dims[] = {n, m, l};
   kiss_fftnd_cfg cfg;

   cfg = kiss_fftnd_alloc(dims, 3, 0, NULL, NULL);
   kiss_fftnd(cfg, in, out);
   free(cfg);
}

void ifft3d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n, int l) {
   int dims[] = {n, m, l};
   kiss_fftnd_cfg cfg;

   cfg = kiss_fftnd_alloc(dims, 3, 1, NULL, NULL);
   kiss_fftnd(cfg, in, out);
   for(int i=0;i<m*n*l;i++) {
      C_MULBYSCALAR(out[i], 1.0/(m*n*l));
   }
   free(cfg);
}
