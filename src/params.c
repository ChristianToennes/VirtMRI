#include "params.h"
#include "stdlib.h"

struct Params* make_params(enum Sequence sequence, int n_params, float* s_params, int xdim, int ydim, int zdim, int xstart, int ystart, int zstart, enum Nearest nearest, bool use_cs, bool use_fft3, struct CSParams* cs_params, struct NoiseParams* noise_params, struct FilterParams* filter_params) {
    struct Params* p = (struct Params*)malloc(sizeof(struct Params));
    p->sequence = sequence;
    p->n_params = n_params;
    p->s_params = s_params;
    p->xdim = xdim;
    p->ydim = ydim;
    p->zdim = zdim;
    p->xstart = xstart;
    p->ystart = ystart;
    p->zstart = zstart;
    p->nearest = nearest;
    p->use_cs = use_cs;
    p->use_fft3 = use_fft3;
    p->cs_params = cs_params;
    p->noise_params = noise_params;
    p->filter_params = filter_params;
    //fprintf(stdout, "%i %i %i %i %i %i %i %i %i %i %i\n", sequence, n_params, s_params, xdim, ydim, zdim, nearest, use_cs, use_fft3, cs_params, noise_params);
    return p;
}

void free_params(struct Params* p) {
    free_noise_params(p->noise_params);
    free_filter_params(p->filter_params);
    free_cs_params(p->cs_params);
    free(p->s_params);
    free(p);
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

void free_dataset(struct Dataset* ds) {
    free(ds->pd);
    free(ds->t1);
    free(ds->t2);
    free(ds->t2s);
    free(ds->na_mm);
    free(ds->na_t1);
    free(ds->na_ex_frac);
    free(ds->na_t2s);
    free(ds->na_t2f);
    free(ds);
}

struct CSParams* make_cs_params(int ninner, int nbreg, double lambda, double lambda2, double mu, double gam, int callback) {
    struct CSParams* params = malloc(sizeof(struct CSParams));
    params->ninner = ninner;
    params->nbreg = nbreg;
    params->lambda = lambda;
    params->lambda2 = lambda2;
    params->mu = mu;
    params->gam = gam;
    params->callback = callback;
    //fprintf(stdout, "%i %i %i %i %i %f %f %f %f %i %f\n", xdim, ydim, zdim, ninner, nbreg, lambda, lambda2, mu, gam, filter_mode, filter_fraction);
    return params;
}

void free_cs_params(struct CSParams* params) {
    free(params);
}

struct NoiseParams* make_noise_params(enum NoiseType noise, double mean, double sigma) {
    struct NoiseParams *params = malloc(sizeof(struct NoiseParams));
    params->noise = noise;
    params->mean = mean;
    params->sigma = sigma;
    //fprintf(stdout, "%i %f %f\n", noise, mean, sigma);
    return params;
}

void free_noise_params(struct NoiseParams* params) {
    free(params);
}

struct FilterParams* make_filter_params(enum FilterMode filter_mode, double filter_fraction, double fmin, double fmax) {
    struct FilterParams* params = malloc(sizeof(struct FilterParams));
    params->filter_mode = filter_mode;
    params->filter_fraction = filter_fraction;
    params->fmin = fmin;
    params->fmax = fmax;
    return params;
}

void free_filter_params(struct FilterParams* params) {
    free(params);
}