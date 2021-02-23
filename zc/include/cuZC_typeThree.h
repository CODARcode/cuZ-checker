#ifndef CUZC_TYPETHREE_H
#define CUZC_TYPETHREE_H

#include "cuZC_typeOne.h"

__global__ void type_three(float *data1, float *data2, double *results, int r3, int r2, int r1, int ssimSize, int ssimShift, int yNum); 
__global__ void gridR_typeThree(double *results, int size); 

#endif /* ----- #ifndef CUZC_TYPETHREE_H  ----- */
