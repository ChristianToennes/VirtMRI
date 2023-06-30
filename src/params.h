#ifndef __params_h
#define __params_h
#include "stdbool.h"
#include <complex.h>

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
void free_noise_params(struct NoiseParams* params);

typedef enum FilterMode {
    Full = 0,
    PseudoRandomLines = 1,
    PseudoRandomPoints = 2,
    RegularLines = 3,
    RandomLines = 4,
    RandomPoints = 5
} filter_mode;

typedef struct FilterParams {
    enum FilterMode filter_mode;
    double filter_fraction;
    double fmin;
    double fmax;
} filter_params;

struct FilterParams* make_filter_params(enum FilterMode filter_mode, double filter_fraction, double fmin, double fmax);
void free_filter_params(struct FilterParams* params);

typedef struct CSParams {
    int ninner;
    int nbreg;
    double lambda;
    double lambda2;
    double mu;
    double gam;
    int callback;
} cs_params;

struct CSParams* make_cs_params(int ninner, int nbreg, double lambda, double lambda2, double mu, double gam, int callback);
void free_cs_params(struct CSParams* params);

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
void free_dataset(struct Dataset* ds);

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
    int ixdim;
    int iydim;
    int izdim;
    int xstart;
    int ystart;
    int zstart;
    enum Nearest nearest;
    bool use_cs;
    bool use_fft3;
    struct CSParams* cs_params;
    struct NoiseParams* noise_params;
    struct FilterParams* filter_params;
} params;

struct Params* make_params(enum Sequence sequence, int n_params, float* s_params, int xdim, int ydim, int zdim, int ixdim, int iydim, int izdim, int xstart, int ystart, int zstart, enum Nearest nearest, bool use_cs, bool use_fft3, struct CSParams* cs_params, struct NoiseParams* noise_params, struct FilterParams* filter_params);
void free_params(struct Params* params);

#endif