#ifndef __tinycs_h
#define __tinycs_h

#include "params.h"

typedef void cs_callback(int z);

void print_stats(char* name, complex float* data, int length);

void compressed_sensing(complex float *f_data, complex float *out, struct Params *params);

#endif