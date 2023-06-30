#ifndef __noise_h
#define __noise_h

#include "params.h"

void addImageNoise(complex float* image, struct Params *p);
bool addKSpaceNoise(complex float* kspace, struct Params *p);

void addGaussianNoise(complex float* image, int img_len, struct NoiseParams *params);

#endif