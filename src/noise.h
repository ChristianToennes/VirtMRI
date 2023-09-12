#ifndef __noise_h
#define __noise_h

#include "params.h"

void addImageNoise(complex float* image, const long dims[], struct Params *p);
bool addKSpaceNoise(complex float* kspace, const long dims[], struct Params *p);

void addGaussianNoise(complex float* image, const long img_len, struct NoiseParams *params);

#endif