#ifndef __tinycs_h
#define __tinycs_h

#include "kissfft/kiss_fft.h"
#include "params.h"

typedef void cs_callback(int z);

void print_stats(char* name, kiss_fft_cpx* data, int length);

void compressed_sensing(kiss_fft_cpx *f_data, kiss_fft_cpx *out, struct Params *params);

#endif