#include "params.h"
#include "stdlib.h"

struct Params* make_params(enum Sequence sequence, int n_params, float* s_params, int xdim, int ydim, int zdim, enum Nearest nearest, bool use_cs, bool use_fft3, struct CSParams* cs_params, struct NoiseParams* noise_params) {
    struct Params* p = (struct Params*)malloc(sizeof(struct Params));
    p->sequence = sequence;
    p->n_params = n_params;
    p->s_params = s_params;
    p->xdim = xdim;
    p->ydim = ydim;
    p->zdim = zdim;
    p->nearest = nearest;
    p->use_cs = use_cs;
    p->use_fft3 = use_fft3;
    p->cs_params = cs_params;
    p->noise_params = noise_params;
    //fprintf(stdout, "%i %i %i %i %i %i %i %i %i %i %i\n", sequence, n_params, s_params, xdim, ydim, zdim, nearest, use_cs, use_fft3, cs_params, noise_params);
    return p;
}

struct Dataset* make_dataset(int kxdim, int kydim, int kzdim, float* pd, float* t1, float* t2, float* t2s, float* na_mm, float* na_t1, float* na_ex_frac, float* na_t2s, float* na_t2f) {
    struct Dataset* ds = (struct Dataset*)malloc(sizeof(struct Dataset));
    ds->k_xdim = kxdim;
    ds->k_ydim = kydim;
    ds->k_zdim = kzdim;
    ds->pd = pd;
    ds->t1 = t1;
    ds->t2 = t2;
    ds->t2s = t2s;
    ds->na_mm = na_mm;
    ds->na_t1 = na_t1;
    ds->na_ex_frac = na_ex_frac;
    ds->na_t2s = na_t2s;
    ds->na_t2f = na_t2f;
    return ds;
}

struct CSParams* make_cs_params(bool filter_only, int xdim, int ydim, int zdim, int ninner, int nbreg, double lambda, double lambda2, double mu, double gam, enum FilterMode filter_mode, double filter_fraction, int callback) {
    struct CSParams* params = malloc(sizeof(struct CSParams));
    params->filter_only = filter_only;
    params->xdim = xdim;
    params->ydim = ydim;
    params->zdim = zdim;
    params->ninner = ninner;
    params->nbreg = nbreg;
    params->lambda = lambda;
    params->lambda2 = lambda2;
    params->mu = mu;
    params->gam = gam;
    params->filter_mode = filter_mode;
    params->filter_fraction = filter_fraction;
    params->callback = callback;
    //fprintf(stdout, "%i %i %i %i %i %f %f %f %f %i %f\n", xdim, ydim, zdim, ninner, nbreg, lambda, lambda2, mu, gam, filter_mode, filter_fraction);
    return params;
}

struct NoiseParams* make_noise_params(enum NoiseType noise, double mean, double sigma) {
    struct NoiseParams *params = malloc(sizeof(struct NoiseParams));
    params->noise = noise;
    params->mean = mean;
    params->sigma = sigma;
    //fprintf(stdout, "%i %f %f\n", noise, mean, sigma);
    return params;
}