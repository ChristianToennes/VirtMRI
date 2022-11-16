#include "kspace_filter.h"

#include "kissfft/_kiss_fft_guts.h"

void apply_pseudo_random_filter(kiss_fft_cpx* kspace, kiss_fft_cpx* out_kspace, struct Params *params) {
    int x,y,z,pos;
    double r,rn,a,b,c;
    bool filter = false;
    int xdim = params->xdim;
    int ydim = params->ydim;
    int zdim = params->zdim;
    int count = 0;
    int count2 = 0;
    double in_frac = 0.1;
    if((params->filter_params->filter_fraction-in_frac) < 0.5) {
        // p = b*(1-in_frac)*0.5
        b = (params->filter_params->filter_fraction-in_frac)/(0.5*(1.0-in_frac));
        a = b/(1.0-in_frac);
    } else {
        // p = (c+b)*0.5*(1-in_frac)
        b = 1.0;
        c = (params->filter_params->filter_fraction-in_frac)/((1.0-in_frac)*0.5) - 1.0;
        a = (b-c)/(1.0-in_frac);
    }
    //fprintf(stderr, "%f %f %f\n", params->filter_fraction, a, b);
    for(z=0;z<zdim;z++) {
        for(y=0;y<ydim;y++) {
            r = (double)y/(double)params->ydim;
            if(r>=0.5) {
                r = 1-r;
            }
            r *= 2;
            if(r<in_frac) {
                filter = false;
            } else {
                rn = (double)rand()/(double)RAND_MAX;
                //filter = pow((1.0-r), params->filter_fraction/0.9) > rn;
                // f(1) = 0
                // f(0.1) = 1
                // a = -1/0.9
                // params->filter_fraction-0.1 = 0.5 * 0.9 * a
                // a = (params->filter_fraction-0.1)/(0.5*0.9)
                // f(x) = -a*x + c
                // 0 = -a*1 + c
                // f(x) = -a*x + a
                filter = rn > (-a*(r-in_frac)+b);
            }
            for(x=0;x<xdim;x++) {
                pos = x+y*xdim+z*xdim*ydim;
                if(filter) {
                    count2++;
                    ASSIGN(out_kspace[pos], 0, 0);
                } else {
                    count++;
                    out_kspace[pos] = kspace[pos];
                }
            }
        }
    }
    //fprintf(stdout, "pseudo randomly sampled: %f discarded: %f\n", (double)count/((double)zdim*(double)ydim*(double)xdim),(double)count2/((double)zdim*(double)ydim*(double)xdim));
}

void apply_pseudo_random_points_filter(kiss_fft_cpx* kspace, kiss_fft_cpx* out_kspace, struct Params *params) {
    int x,y,z,pos;
    double r,rn,a,b,c;
    bool filter = false;
    int xdim = params->xdim;
    int ydim = params->ydim;
    int zdim = params->zdim;
    int count = 0;
    int count2 = 0;
    double in_frac = 0.1;
    if((params->filter_params->filter_fraction-in_frac) < 0.5) {
        // p = b*(1-in_frac)*0.5
        b = (params->filter_params->filter_fraction-in_frac)/(0.5*(1.0-in_frac));
        a = b/(1.0-in_frac);
    } else {
        // p = (c+b)*0.5*(1-in_frac)
        b = 1.0;
        c = (params->filter_params->filter_fraction-in_frac)/((1.0-in_frac)*0.5) - 1.0;
        a = (b-c)/(1.0-in_frac);
    }
    //fprintf(stderr, "%f %f %f\n", params->filter_fraction, a, b);
    for(z=0;z<zdim;z++) {
        for(y=0;y<ydim;y++) {
            for(x=0;x<xdim;x++) {
                if(y > 0.5*(double)params->ydim) {
                    if(x > 0.5*(double)params->xdim) {
                        r = sqrt(((double)params->ydim-(double)y)*((double)params->ydim-(double)y)+((double)params->xdim-(double)x)*((double)params->xdim-(double)x))/(double)params->ydim;
                    } else {
                        r = sqrt(((double)params->ydim-(double)y)*((double)params->ydim-(double)y)+(double)x*(double)x)/(double)params->ydim;
                    }
                } else {
                    if(x > 0.5*(double)params->xdim) {
                        r = sqrt((double)y*(double)y+((double)params->xdim-(double)x)*((double)params->xdim-(double)x))/(double)params->ydim;
                    } else {
                        r = sqrt((double)y*(double)y+(double)x*x)/(double)params->ydim;
                    }
                }
                r *= 2;
                if(r<in_frac) {
                    filter = false;
                } else {
                    rn = (double)rand()/(double)RAND_MAX;
                    //filter = pow((1.0-r), params->filter_fraction/0.9) > rn;
                    // f(1) = 0
                    // f(0.1) = 1
                    // a = -1/0.9
                    // params->filter_fraction-0.1 = 0.5 * 0.9 * a
                    // a = (params->filter_fraction-0.1)/(0.5*0.9)
                    // f(x) = -a*x + c
                    // 0 = -a*1 + c
                    // f(x) = -a*x + a
                    filter = rn > (-a*(r-in_frac)+b);
                }
                pos = x+y*xdim+z*xdim*ydim;
                if(filter) {
                    count2++;
                    ASSIGN(out_kspace[pos], 0, 0);
                } else {
                    count++;
                    out_kspace[pos] = kspace[pos];
                }
            }
        }
    }
    //fprintf(stdout, "pseudo randomly sampled: %f discarded: %f\n", (double)count/((double)zdim*(double)ydim*(double)xdim),(double)count2/((double)zdim*(double)ydim*(double)xdim));
}

void apply_random_filter(kiss_fft_cpx* kspace, kiss_fft_cpx* out_kspace, struct Params *params) {
    int count = 0;
    int count2 = 0;
    int x,y,z,pos;
    int xdim = params->xdim;
    int ydim = params->ydim;
    int zdim = params->zdim;
    double rn;
    for(z=0;z<zdim;z++) {
        for(y=0;y<ydim;y++) {
            for(x=0;x<xdim;x++) {    
                rn = (double)rand()/(double)RAND_MAX;
                pos = x+y*xdim+z*xdim*ydim;
                if(rn < params->filter_params->filter_fraction == 0) {
                    count2++;
                    ASSIGN(out_kspace[pos], 0, 0);
                } else {
                    count++;
                    out_kspace[pos] = kspace[pos];
                }
            }
        }
    }
    //fprintf(stdout, "randomly sampled: %f discarded: %f\n", (double)count/((double)zdim*(double)ydim*(double)xdim),(double)count2/((double)zdim*(double)ydim*(double)xdim));
}

void apply_random_lines_filter(kiss_fft_cpx* kspace, kiss_fft_cpx* out_kspace, struct Params *params) {
    int count = 0;
    int count2 = 0;
    int x,y,z,pos;
    int xdim = params->xdim;
    int ydim = params->ydim;
    int zdim = params->zdim;
    double rn;
    for(z=0;z<zdim;z++) {
        for(y=0;y<ydim;y++) {
            rn = (double)rand()/(double)RAND_MAX;
            for(x=0;x<xdim;x++) {    
                pos = x+y*xdim+z*xdim*ydim;
                if(rn < params->filter_params->filter_fraction == 0) {
                    count2++;
                    ASSIGN(out_kspace[pos], 0, 0);
                } else {
                    count++;
                    out_kspace[pos] = kspace[pos];
                }
            }
        }
    }
    //fprintf(stdout, "random lines sampled: %f discarded: %f\n", (double)count/((double)zdim*(double)ydim*(double)xdim),(double)count2/((double)zdim*(double)ydim*(double)xdim));
}

void apply_regular_filter(kiss_fft_cpx* kspace, kiss_fft_cpx* out_kspace, struct Params *params) {
    int count = 0;
    int count2 = 0;
    int x,y,z,pos;
    int xdim = params->xdim;
    int ydim = params->ydim;
    int zdim = params->zdim;
    bool filtered = true;
    double nth = params->filter_params->filter_fraction;
    if(params->filter_params->filter_fraction > 0.5) {
        nth = 1.0-params->filter_params->filter_fraction;
    }
    for(z=0;z<zdim;z++) {
        for(y=0;y<ydim;y++) {
            filtered = fabs(round(y*nth)/nth-y) < 1;
            if(params->filter_params->filter_fraction > 0.5) {
                filtered = !filtered;
            }
            //fprintf(stdout, "%i %f %f\n", y, round(y*100/nth)*nth/100.0, fabs(round(y*100/nth)*nth/100.0-y));
            for(x=0;x<xdim;x++) {
                pos = x+y*xdim+z*xdim*ydim;
                if(filtered) {
                    count2++;
                    ASSIGN(out_kspace[pos], 0, 0);
                } else {
                    count++;
                    out_kspace[pos] = kspace[pos];
                }
            }
        }
    }
    //fprintf(stdout, "regularly spaced sampled: %f discarded: %f frac: %f nth: %f\n", (double)count/((double)zdim*(double)ydim*(double)xdim),(double)count2/((double)zdim*(double)ydim*(double)xdim), params->filter_fraction, nth);
}

void apply_freq_cutoff(kiss_fft_cpx* kspace, kiss_fft_cpx* out_kspace, struct Params *params) {
    if(params->use_fft3) {
        fprintf(stdout, "%f %f\n", params->filter_params->fmin, params->filter_params->fmax);
        double maxf = sqrt(0.5*params->xdim*0.5*params->xdim+0.5*params->ydim*0.5*params->ydim+0.5*params->zdim*0.5*params->zdim);
        double fmin = params->filter_params->fmin*0.01*maxf;
        double fmax = params->filter_params->fmax*0.01*maxf;
        double d = 0;
        fprintf(stdout, "%f %f %f\n", fmin, fmax, maxf);
        for(int x=0;x<params->xdim/2;x++) {
            
            for(int y=0;y<params->ydim/2;y++) {
                for(int z=0;z<params->zdim/2;z++) {
                    d = sqrt(x*x+y*y+z*z);
                    int i = x+y*params->xdim+z*params->ydim*params->xdim;
                    if(i>=256*256*256) fprintf(stdout, "1 %i %i %i %i %i\n", x,y,z,i,256*256*256);
                    if(d<fmin || d>fmax) {
                        ASSIGN(out_kspace[i], 0, 0);
                        if(i>=256*256*256) fprintf(stdout, "2 %i %i %i %i %i\n", x,y,z,i,256*256*256);
                        i = (params->xdim-x-1)+y*params->xdim+z*params->ydim*params->xdim;
                        ASSIGN(out_kspace[i], 0, 0);
                        if(i>=256*256*256) fprintf(stdout, "3 %i %i %i %i %i\n", x,y,z,i,256*256*256);
                        i = x+(params->ydim-y-1)*params->xdim+z*params->ydim*params->xdim;
                        ASSIGN(out_kspace[i], 0, 0);
                        if(i>=256*256*256) fprintf(stdout, "4 %i %i %i %i %i\n", x,y,z,i,256*256*256);
                        i = x+y*params->xdim+(params->zdim-z-1)*params->ydim*params->xdim;
                        ASSIGN(out_kspace[i], 0, 0);
                        if(i>=256*256*256) fprintf(stdout, "5 %i %i %i %i %i\n", x,y,z,i,256*256*256);

                        i = (params->xdim-x-1)+(params->ydim-y-1)*params->xdim+z*params->ydim*params->xdim;
                        ASSIGN(out_kspace[i], 0, 0);
                        if(i>=256*256*256) fprintf(stdout, "6 %i %i %i %i %i\n", x,y,z,i,256*256*256);
                        i = x+(params->ydim-y-1)*params->xdim+(params->zdim-z-1)*params->ydim*params->xdim;
                        ASSIGN(out_kspace[i], 0, 0);
                        if(i>=256*256*256) fprintf(stdout, "7 %i %i %i %i %i\n", x,y,z,i,256*256*256);
                        i = (params->xdim-x-1)+y*params->xdim+(params->zdim-z-1)*params->ydim*params->xdim;
                        ASSIGN(out_kspace[i], 0, 0);
                        if(i>=256*256*256) fprintf(stdout, "8 %i %i %i %i %i\n", x,y,z,i,256*256*256);

                        i = (params->xdim-x-1)+(params->ydim-y-1)*params->xdim+(params->zdim-z-1)*params->ydim*params->xdim;
                        ASSIGN(out_kspace[i], 0, 0);
                        if(i>=256*256*256) fprintf(stdout, "9 %i %i %i %i %i\n", x,y,z,i,256*256*256);
                        
                    } else {
                        out_kspace[i] = kspace[i];
                        if(i>=256*256*256) fprintf(stdout, "21 %i %i %i %i %i\n", x,y,z,i,256*256*256);
                        i = (params->xdim-x-1)+y*params->xdim+z*params->ydim*params->xdim;
                        out_kspace[i] = kspace[i];
                        if(i>=256*256*256) fprintf(stdout, "22 %i %i %i %i %i\n", x,y,z,i,256*256*256);
                        i = x+(params->ydim-y-1)*params->xdim+z*params->ydim*params->xdim;
                        out_kspace[i] = kspace[i];
                        if(i>=256*256*256) fprintf(stdout, "23 %i %i %i %i %i\n", x,y,z,i,256*256*256);
                        i = x+y*params->xdim+(params->zdim-z-1)*params->ydim*params->xdim;
                        out_kspace[i] = kspace[i];
                        if(i>=256*256*256) fprintf(stdout, "24 %i %i %i %i %i\n", x,y,z,i,256*256*256);

                        i = (params->xdim-x-1)+(params->ydim-y-1)*params->xdim+z*params->ydim*params->xdim;
                        out_kspace[i] = kspace[i];
                        if(i>=256*256*256) fprintf(stdout, "25 %i %i %i %i %i\n", x,y,z,i,256*256*256);
                        i = x+(params->ydim-y-1)*params->xdim+(params->zdim-z-1)*params->ydim*params->xdim;
                        out_kspace[i] = kspace[i];
                        if(i>=256*256*256) fprintf(stdout, "26 %i %i %i %i %i\n", x,y,z,i,256*256*256);
                        i = (params->xdim-x-1)+y*params->xdim+(params->zdim-z-1)*params->ydim*params->xdim;
                        out_kspace[i] = kspace[i];
                        if(i>=256*256*256) fprintf(stdout, "27 %i %i %i %i %i\n", x,y,z,i,256*256*256);

                        i = (params->xdim-x-1)+(params->ydim-y-1)*params->xdim+(params->zdim-z-1)*params->ydim*params->xdim;
                        out_kspace[i] = kspace[i];
                        if(i>=256*256*256) fprintf(stdout, "28 %i %i %i %i %i\n", x,y,z,i,256*256*256);
                    }
                }
            }
        }
    }
    else {
        double maxf = params->xdim*params->xdim+params->ydim*params->ydim;
        double fmin = params->filter_params->fmin*0.01*maxf;
        double fmax = params->filter_params->fmax*0.01*maxf;
        for(int x=0;x<params->xdim/2;x++) {
            for(int y=0;y<params->ydim/2;y++) {
                double d = x*x+y*y;
                for(int z=0;z<params->zdim;z++) {
                    int i = x+y*params->xdim+z*params->ydim*params->xdim;
                    if(d<fmin || d>fmax) {
                        ASSIGN(out_kspace[i], 0, 0);
                        i = (params->xdim-x-1)+y*params->xdim+z*params->ydim*params->xdim;
                        ASSIGN(out_kspace[i], 0, 0);
                        i = x+(params->ydim-y-1)*params->xdim+z*params->ydim*params->xdim;
                        ASSIGN(out_kspace[i], 0, 0);

                        i = (params->xdim-x-1)+(params->ydim-y-1)*params->xdim+z*params->ydim*params->xdim;
                        ASSIGN(out_kspace[i], 0, 0);
                    } else {
                        out_kspace[i] = kspace[i];
                        i = params->xdim-x+y*params->xdim+z*params->ydim*params->xdim;
                        out_kspace[i] = kspace[i];
                        i = x+(params->ydim-y)*params->xdim+z*params->ydim*params->xdim;
                        out_kspace[i] = kspace[i];

                        i = params->xdim-x+(params->ydim-y)*params->xdim+z*params->ydim*params->xdim;
                        out_kspace[i] = kspace[i];
                    }
                }
            }
        }
    }
}

bool apply_kspace_filter(kiss_fft_cpx* kspace, kiss_fft_cpx* out_kspace, struct Params *params) {
    bool modified = false;
    switch(params->filter_params->filter_mode) {
        case PseudoRandomLines:
            modified = true;
            apply_pseudo_random_filter(kspace, out_kspace, params);
            break;
        case RegularLines:
            modified = true;
            apply_regular_filter(kspace, out_kspace, params);
            break;
        case RandomPoints:
            modified = true;
            apply_random_filter(kspace, out_kspace, params);
            break;
        case RandomLines:
            modified = true;
            apply_random_lines_filter(kspace, out_kspace, params);
            break;
        case PseudoRandomPoints:
            modified = true;
            apply_pseudo_random_points_filter(kspace, out_kspace, params);
            break;
        default:
            modified = false;
            break;
    }
    if(params->filter_params->fmin > 0 || params->filter_params->fmax < 100) {
        modified = true;
        apply_freq_cutoff(kspace, out_kspace, params);
    }

    return modified;
}
