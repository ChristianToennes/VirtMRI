#ifndef __noise_h
#define __noise_h

#include "kissfft/kiss_fft.h"
#include "params.h"

void addImageNoise(kiss_fft_cpx* image, struct Params *p);
bool addKSpaceNoise(kiss_fft_cpx* kspace, struct Params *p);

void addGaussianNoise(kiss_fft_cpx* image, int img_len, struct NoiseParams *params);

#endif