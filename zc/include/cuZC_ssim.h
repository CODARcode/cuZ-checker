#ifndef CUZC_SSIM_H
#define CUZC_SSIM_H

#include <stdio.h>
 
const int N = 16; 
const int blocksize = 16; 

#define FULL_MASK 0xffffffff
 
__inline__ __device__
int warpReduceSum(int val) {

    for (int offset = warpSize/2; offset > 0; offset /= 2) 
        //val += __shfl_down(val, offset);
        val += __shfl_down_sync(FULL_MASK, val, offset);
        //val += 1;
    return val;
  
}

__inline__ __device__
int blockReduceSum(int val) {

    static __shared__ int shared[32]; // Shared mem for 32 partial sums
    int lane = threadIdx.x % warpSize;
    int wid = threadIdx.x / warpSize;
    printf("warp=%i\n", blockDim.x);

    val = warpReduceSum(val);     // Each warp performs partial reduction

    if (lane==0) shared[wid]=val; // Write reduced value to shared memory

    __syncthreads();              // Wait for all partial reductions

    //read from shared memory only if that warp existed
    val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;

    if (wid==0) val = warpReduceSum(val); //Final reduce within first warp

    return val;
    
}

__device__
double ssim_winComp(float xMin, float xMax, float yMin, float yMax, float xSum, float x2Sum, float ySum, float y2Sum, float xySum, double np);

__global__ void hello(char *a, int *b); 

__global__ void ssim(float *data1, float *data2, double *results, int r3, int r2, int r1, int ssimSize, int ssimShift); 

#endif /* ----- #ifndef CUZC_SSIM_H  ----- */
