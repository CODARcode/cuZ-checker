#ifndef CUZC_TYPEONE_H
#define CUZC_TYPEONE_H

#include <cuda_runtime.h>

// Utilities and system includes
#include <helper_cuda.h>  // helper function CUDA error checking and initialization
#include <helper_functions.h>  // helper for shared functions common to CUDA Samples

__device__
void reduction(double sum1, double sum2,
        double minDiff, double maxDiff, double sumDiff, double sumOfDiffSquare, 
        double minErr, double maxErr, double sumErr, double sumErrSqr);

__global__ void type_one(float *data1, float *data2, double *diff, double *results, int r3, int r2, int r1, size_t ne);

__global__ void gridReduction(double *results, int r3); 

#endif /* ----- #ifndef CUZC_TYPEONE_H  ----- */
