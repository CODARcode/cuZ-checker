#ifndef CUZC_TYPETWO_H
#define CUZC_TYPETWO_H

#include "cuZC_typeOne.h"

__global__ void type_two(float *data, float *der, float *autocor, int r3, int r2, int r1, float avg, size_t order); 

#endif /* ----- #ifndef CUZC_TYPETWO_H  ----- */
