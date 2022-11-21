#include "coords.h"

int static inline toImgX(int x, struct Params* p) {
    return p->xstart + x;
}

int static inline toImgY(int y, struct Params* p) {
    return p->ystart + y;
}

int static inline toImgZ(int z, struct Params* p) {
    return p->zstart + z;
}

int static inline toDSX(int x, struct Params* p) {
    return x - p->xstart;
}

int static inline toDSY(int y, struct Params* p) {
    return y - p->ystart;
}

int static inline toDSZ(int z, struct Params* p) {
    return z - p->zstart;
}
