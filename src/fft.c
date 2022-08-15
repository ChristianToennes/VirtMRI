#include "fft.h"

#include <stdint.h>
#include <stdio.h>

#include "_kiss_fft_guts.h"

void fft(kiss_fft_cpx *in, kiss_fft_cpx *out, int m) {
   kiss_fft_cfg cfg;

   cfg = kiss_fft_alloc(m, 0, NULL, NULL);
   kiss_fft(cfg, in, out);
   free(cfg);
}

void ifft(kiss_fft_cpx *in, kiss_fft_cpx *out, int m) {
   kiss_fft_cfg cfg;

   cfg = kiss_fft_alloc(m, 1, NULL, NULL);
   kiss_fft(cfg, in, out);
   for(int i=0;i<m;i++) {
      C_MULBYSCALAR(out[i], 1.0/m);
   }
   free(cfg);
}

void fft_stride(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int nl) {
   for(int i=0;i<nl;i++) {
      kiss_fft_cfg cfg;
      cfg = kiss_fft_alloc (m, 0, NULL, NULL);
      kiss_fft_stride(cfg, in+i, out+i*nl, nl);
      free(cfg);
   }
}