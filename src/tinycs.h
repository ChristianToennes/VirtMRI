#ifndef __tinycs_h
#define __tinycs_h

#include "kissfft/kiss_fft.h"
#include "params.h"

typedef void cs_callback(int z);

void print_stats(char* name, kiss_fft_cpx* data, int length);

void apply_pseudo_random_filter(kiss_fft_cpx* kspace, kiss_fft_cpx* out_kspace, struct CSParams *params);

void apply_regular_filter(kiss_fft_cpx* kspace, kiss_fft_cpx* out_kspace, struct CSParams *params);

void apply_random_filter(kiss_fft_cpx* kspace, kiss_fft_cpx* out_kspace, struct CSParams *params);

void apply_cs_filter(kiss_fft_cpx* kspace, kiss_fft_cpx* out_kspace, struct CSParams *params);

void compressed_sensing(kiss_fft_cpx *f_data, kiss_fft_cpx *out, struct CSParams *params);

#endif