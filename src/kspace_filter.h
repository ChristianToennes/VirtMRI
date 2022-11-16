#ifndef __kspace_filter_h
#define __kspace_filter_h

#include "kissfft/kiss_fft.h"
#include "params.h"

void apply_pseudo_random_filter(kiss_fft_cpx* kspace, kiss_fft_cpx* out_kspace, struct Params *params);

void apply_regular_filter(kiss_fft_cpx* kspace, kiss_fft_cpx* out_kspace, struct Params *params);

void apply_random_filter(kiss_fft_cpx* kspace, kiss_fft_cpx* out_kspace, struct Params *params);

void apply_freq_cutoff(kiss_fft_cpx* kspace, kiss_fft_cpx* out_kspace, struct Params *params);

bool apply_kspace_filter(kiss_fft_cpx* kspace, kiss_fft_cpx* out_kspace, struct Params *params);

#endif