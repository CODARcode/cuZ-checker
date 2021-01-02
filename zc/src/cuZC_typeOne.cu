#include <stdio.h>
#include <math.h>
#include "cuZC_ssim.h"
#include "cuZC_typeOne.h"
#include "matrix.hpp"

__device__
void reduction(double sum1, double sum2,
        double minDiff, double maxDiff, double sumDiff, double sumOfDiffSquare, 
        double minErr, double maxErr, double sumErr, double sumErrSqr, double *results){

    static __shared__ double shared[32*10];

    int lane = threadIdx.x;
    int wid = threadIdx.y;


    for (int offset = warpSize/2; offset > 0; offset /= 2) 
    {
        minDiff = min(minDiff, __shfl_xor_sync(FULL_MASK, minDiff, offset));
        maxDiff = max(maxDiff, __shfl_xor_sync(FULL_MASK, maxDiff, offset));
        minErr = min(minErr, __shfl_xor_sync(FULL_MASK, minErr, offset));
        maxErr = max(maxErr, __shfl_xor_sync(FULL_MASK, maxErr, offset));
        sum1 += __shfl_down_sync(FULL_MASK, sum1, offset);
        sum2 += __shfl_down_sync(FULL_MASK, sum2, offset);
        sumDiff += __shfl_down_sync(FULL_MASK, sumDiff, offset);
        sumOfDiffSquare += __shfl_down_sync(FULL_MASK, sumOfDiffSquare, offset);
        sumErr += __shfl_down_sync(FULL_MASK, sumErr, offset);
        sumErrSqr += __shfl_down_sync(FULL_MASK, sumErrSqr, offset);
    }

    if (lane==0){
        shared[wid] = minDiff;
        shared[32+wid] = maxDiff;
        shared[32*2+wid] = minErr;
        shared[32*3+wid] = maxErr;
        shared[32*4+wid] = sum1;
        shared[32*5+wid] = sum2;
        shared[32*6+wid] = sumDiff;
        shared[32*7+wid] = sumOfDiffSquare;
        shared[32*8+wid] = sumErr;
        shared[32*9+wid] = sumErrSqr;
    }

    __syncthreads();                  

    //if (wid==0)printf("ddata%i=%e:%e\n", 32*6+lane, shared[32*6+lane], ySum);

    if (wid==0){
        if (threadIdx.x < blockDim.y){
            minDiff = shared[lane];
            maxDiff = shared[32+lane];
            minErr = shared[32*2+lane];
            maxErr = shared[32*3+lane];
            sum1 = shared[32*4+lane];
            sum2 = shared[32*5+lane];
            sumDiff = shared[32*6+lane];
            sumOfDiffSquare = shared[32*7+lane];
            sumErr = shared[32*8+lane];
            sumErrSqr = shared[32*9+lane];
        }else{
            minDiff = shared[0];  
            maxDiff = shared[32]; 
            minErr = shared[32*2]; 
            maxErr = shared[32*3]; 
            sum1 = 0; 
            sum2 = 0;
            sumDiff = 0; 
            sumOfDiffSquare = 0;
            sumErr = 0;
            sumErrSqr = 0;
        }

        for (int offset = warpSize/2; offset > 0; offset /= 2) 
        {
            minDiff = min(minDiff, __shfl_xor_sync(FULL_MASK, minDiff, offset));
            maxDiff = max(maxDiff, __shfl_xor_sync(FULL_MASK, maxDiff, offset));
            minErr = min(minErr, __shfl_xor_sync(FULL_MASK, minErr, offset));
            maxErr = max(maxErr, __shfl_xor_sync(FULL_MASK, maxErr, offset));
            sum1 += __shfl_down_sync(FULL_MASK, sum1, offset);
            sum2 += __shfl_down_sync(FULL_MASK, sum2, offset);
            sumDiff += __shfl_down_sync(FULL_MASK, sumDiff, offset);
            sumOfDiffSquare += __shfl_down_sync(FULL_MASK, sumOfDiffSquare, offset);
            sumErr += __shfl_down_sync(FULL_MASK, sumErr, offset);
            sumErrSqr += __shfl_down_sync(FULL_MASK, sumErrSqr, offset);
        }
        
        if (lane==0){
            results[blockIdx.x] = minDiff;
            results[gridDim.x+blockIdx.x] = minErr;
            results[gridDim.x*2+blockIdx.x] = maxDiff;
            results[gridDim.x*3+blockIdx.x] = maxErr;
            results[gridDim.x*4+blockIdx.x] = sum1;
            results[gridDim.x*5+blockIdx.x] = sum2;
            results[gridDim.x*6+blockIdx.x] = sumDiff;
            results[gridDim.x*7+blockIdx.x] = sumOfDiffSquare;
            results[gridDim.x*8+blockIdx.x] = sumErr;
            results[gridDim.x*9+blockIdx.x] = sumErrSqr;
        }
    }
    //if(lane==0){
    //    if (sum1>0.0)printf("test%i,%i,%i:%e\n",lane,wid,blockIdx.x, sum1);

    //}

}

__global__ void type_one(float *data1, float *data2, double *diff, double *results, int r3, int r2, int r1, size_t ne) 
{
    int tidx = threadIdx.x;
    int tidy = threadIdx.y;
    int bid = blockIdx.x;

    float Data1, Data2;
    double Diff;
	double minDiff = data2[0]-data1[0];
	double maxDiff = minDiff;
	double sum1 = 0, sum2 = 0, sumDiff = 0, sumOfDiffSquare = 0; 
	
	double err;
	size_t numOfElem = ne;
	double minErr = fabs(minDiff);
	double maxErr = minErr;
	double sumErr = 0, sumErrSqr = 0;

    int i, j;

    for (j=tidy; j<r2; j+=blockDim.y){
        for (i=tidx; i<r1; i+=blockDim.x){
            Data1 = data1[bid*r1*r2+j*r1+i];
            Data2 = data2[bid*r1*r2+j*r1+i];
            sum1 += Data1; 
            sum2 += Data2;

            Diff = Data2 - Data1;
            diff[bid*r1*r2+j*r1+i] = Diff;
            minDiff = min(minDiff, Diff);
            maxDiff = max(maxDiff, Diff);
            sumDiff += Diff;
            sumOfDiffSquare += Diff * Diff;

            err = fabs(Diff);
            minErr = min(minErr, err);
            maxErr = max(maxErr, err);
            sumErr += err;
            sumErrSqr += err*err;
        }
    }
    __syncthreads();                  

    reduction(sum1, sum2, minDiff, maxDiff, sumDiff, sumOfDiffSquare, minErr, maxErr, sumErr, sumErrSqr, results);

//if (tid == 0)printf("ydata%i,%i=%e\n",Offsetx,Offsety, result);
//    if (tid==0) results[bid] = result;
    
}

__global__ void gridReduction(double *results, int r3) 
{
    int tidx = threadIdx.x;
    int tidy = threadIdx.y;
    double data = results[tidy*r3+tidx];

    for (int i=(tidx+blockDim.x); i<r3; i+=blockDim.x){
        if (tidy<2) data = min(data, results[tidy*r3+i]);
        else if (tidy<4) data = max(data, results[tidy*r3+i]);
        else data += results[tidy*r3+i];
    }

    for (int offset = warpSize/2; offset > 0; offset /= 2) 
    {
        if (tidy<2) data = min(data, __shfl_xor_sync(FULL_MASK, data, offset));
        else if (tidy<4) data = max(data, __shfl_xor_sync(FULL_MASK, data, offset));
        else data += __shfl_down_sync(FULL_MASK, data, offset);
    }

    if (tidx==0) results[tidy*r3] = data;
        
}
