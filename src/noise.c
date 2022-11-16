#include "noise.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include "_kiss_fft_guts.h"
#include "stdbool.h"

bool addKSpaceNoise(kiss_fft_cpx* kspace, struct Params *p) {
    bool modified = false;
    if(p->noise_params->noise & Gaussian) {
        modified = true;
        addGaussianNoise(kspace, p->xdim*p->ydim*p->zdim, p->noise_params);
    }
    return modified;
}

void addImageNoise(kiss_fft_cpx* image, struct Params *p) {

}

void addGaussianNoise(kiss_fft_cpx* image, int image_len, struct NoiseParams *params) {
    double u,v,s,x,y;
    kiss_fft_cpx r;
    for(int i=0;i<image_len;i++) {
        do {
            u = (double)rand()/(double)RAND_MAX;
        } while(u == 0);
        do {
            v = (double)rand()/(double)RAND_MAX;
        } while(v == 0);
        s = sqrt(-2.0*log(u));
        x = s * cos(2.0*M_PI*v);
        y = s * sin(2.0*M_PI*v);
        ASSIGN(r, (params->sigma*x+params->mean), (params->sigma*y+params->mean));
        C_ADDTO(image[i], r);
    }
}
