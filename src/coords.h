#ifndef __coords_h
#define __coords_h

#include "params.h"

int static inline toImgX(int x, struct Params* params);
int static inline toImgY(int y, struct Params* params);
int static inline toImgZ(int z, struct Params* params);

int static inline toDSX(int x, struct Params* params);
int static inline toDSY(int y, struct Params* params);
int static inline toDSZ(int z, struct Params* params);

#endif