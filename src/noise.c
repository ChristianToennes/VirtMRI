#include "noise.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include "stdbool.h"
#include "stdlib.h"

bool addKSpaceNoise(complex float* kspace, const long dims[], struct Params *p) {
    bool modified = false;
    if(p->noise_params->noise & Gaussian) {
        modified = true;
        addGaussianNoise(kspace, dims[0]*dims[1]*dims[2]*dims[3], p->noise_params);
    }
    return modified;
}

void addImageNoise(complex float* image, const long dims[], struct Params *p) {

}

void addGaussianNoise(complex float* image, const long image_len, struct NoiseParams *params) {
    double u,v,s,x,y;
    complex float r;
    for(long i=0;i<image_len;i++) {
        do {
            u = (double)rand()/(double)RAND_MAX;
        } while(u == 0);
        do {
            v = (double)rand()/(double)RAND_MAX;
        } while(v == 0);
        s = sqrt(-2.0*log(u));
        x = s * cos(2.0*M_PI*v);
        y = s * sin(2.0*M_PI*v);
        r = params->sigma*x+params->mean+(params->sigma*y+params->mean)*I;
        image[i] += r;
    }
}
