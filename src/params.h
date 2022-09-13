#ifndef __params_h
#define __params_h
#include "stdbool.h"
#include "kissfft/kiss_fft.h"

void fft2d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n);
void ifft2d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n);
void fft3d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n, int l);
void ifft3d(kiss_fft_cpx *in, kiss_fft_cpx *out, int m, int n, int l);

typedef enum NoiseType {
    Gaussian = 1,
    Motion = 2,
} noise_type;

typedef struct NoiseParams {
    enum NoiseType noise;
    double mean;
    double sigma;
} noise_params;

struct NoiseParams* make_noise_params(enum NoiseType noise, double mean, double sigma);

typedef struct CSParams {
    int xdim;
    int ydim;
    int zdim;
    int ninner;
    int nbreg;
    double lambda;
    double lambda2;
    double mu;
    double gam;
    double filter_fraction;
    int callback;
} cs_params;

struct CSParams* make_cs_params(int xdim, int ydim, int zdim, int ninner, int nbreg, double lambda, double lambda2, double mu, double gam, double filter_fraction, int callback);

typedef struct Dataset {
    float* pd;
    float* t1;
    float* t2;
    float* t2s;
    float* na_mm;
    float* na_t1;
    float* na_ex_frac;
    float* na_t2s;
    float* na_t2f;
    int k_xdim;
    int k_ydim;
    int k_zdim;
} dataset;

struct Dataset* make_dataset(int kxdim, int kydim, int kzdim, float* pd, float* t1, float* t2, float* t2s, float* na_mm, float* na_t1, float* na_ex_frac, float* na_t2s, float* na_t2f);

typedef enum Sequence {
    SE, IR, bSSFP, FISP, PSIF, SGRE, Na, NaSQ, NaTQ, NaTQSQ, NaTQF, pcbSSFP
} sequence;

typedef enum Nearest {
    Nearest, KNearest, Average
} nearest;

typedef struct Params {
    enum Sequence sequence;
    int n_params;
    float* s_params;
    int xdim;
    int ydim;
    int zdim;
    enum Nearest nearest;
    bool use_cs;
    bool use_fft3;
    struct CSParams* cs_params;
    struct NoiseParams* noise_params;
} params;

struct Params* make_params(enum Sequence sequence, int n_params, float* s_params, int xdim, int ydim, int zdim, enum Nearest nearest, bool use_cs, bool use_fft3, struct CSParams* cs_params, struct NoiseParams* noise_params);

#endif