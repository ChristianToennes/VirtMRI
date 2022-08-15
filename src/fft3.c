#include "fft3.h"

#include <stdint.h>
#include <stdio.h>

#include "_kiss_fft_guts.h"

void __fft3d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n, int l) {
   int dims[] = {n, m, l};
   kiss_fftnd_cfg cfg;

   cfg = kiss_fftnd_alloc(dims, 3, 0, NULL, NULL);
   kiss_fftnd(cfg, in, out);
   free(cfg);
}

void __ifft3d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n, int l) {
   int dims[] = {n, m, l};
   kiss_fftnd_cfg cfg;

   cfg = kiss_fftnd_alloc(dims, 3, 1, NULL, NULL);
   kiss_fftnd(cfg, in, out);
   free(cfg);
}


void fft3d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n, int l) {
   int dims[] = {n, m};
   kiss_fftnd_cfg cfg;
   kiss_fft_cfg cfg_stride;

   kiss_fft_cpx min, max, sum;
   double avg;
   min = in[0];
   max = in[0];
   ASSIGN(sum, 0, 0);
   avg=0;
   for(int i=0;i<m*n*l;i++) {
      if (C_GT(min, in[i])) { min = in[i];}
      if (C_LT(max, in[i])) { max = in[i];}
      C_ADDTO(sum, in[i]);
   }
   avg = (REAL(sum)*REAL(sum)+IMAG(sum)*IMAG(sum)) / (double)(m*n*l);
   fprintf(stdout, "in %f %fi %f %fi %f\n", REAL(min), IMAG(min), REAL(max), IMAG(max), avg);

   for(int z=0;z<l;z++) {
      cfg = kiss_fftnd_alloc(dims, 2, 0, NULL, NULL);
      kiss_fftnd(cfg, &in[z*m*n], &out[z*m*n]);
      free(cfg);
   }

   for(int i=0;i<m*n;i++) {
      cfg_stride = kiss_fft_alloc(l, 0, NULL, NULL);
      kiss_fft_stride(cfg_stride, &in[i], &out[i], m*n);
      free(cfg_stride);
   }
}

void ifft3d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n, int l) {
   int dims[] = {n, m};
   kiss_fftnd_cfg cfg;
   kiss_fft_cfg cfg_stride;

   for(int z=0;z<l;z++) {
      cfg = kiss_fftnd_alloc(dims, 2, 1, NULL, NULL);
      kiss_fftnd(cfg, &in[z*m*n], &out[z*m*n]);
      free(cfg);
   }

   for(int i=0;i<m*n;i++) {
      cfg_stride = kiss_fft_alloc(l, 1, NULL, NULL);
      kiss_fft_stride(cfg_stride, &in[i], &out[i], m*n);
      free(cfg_stride);
   }
   for(int i=0;i<m*n;i++) {
      C_MULBYSCALAR(out[i], 1.0/(m*n*l));
   }
}