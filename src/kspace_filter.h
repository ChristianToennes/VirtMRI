#ifndef __kspace_filter_h
#define __kspace_filter_h

#include "params.h"

void apply_pseudo_random_filter(complex float* kspace, complex float* out_kspace, struct Params *params);

void apply_regular_filter(complex float* kspace, complex float* out_kspace, struct Params *params);

void apply_random_filter(complex float* kspace, complex float* out_kspace, struct Params *params);

void apply_freq_cutoff(complex float* kspace, complex float* out_kspace, struct Params *params);

bool apply_kspace_filter(complex float* kspace, complex float* out_kspace, struct Params *params);

#endif