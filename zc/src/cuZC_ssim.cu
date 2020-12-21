#include <stdio.h>
#include <math.h>
#include "cuZC_ssim.h"
#include "matrix.hpp"

__global__ void hello(char *a, int *b) 
{
    a[threadIdx.x] += b[threadIdx.x];
}

__device__
double ssim_winComp(float xMin, float xMax, float yMin, float yMax, float xSum, float x2Sum, float ySum, float y2Sum, float xySum, double np) {

    static __shared__ float shared[32*9]; // Shared mem for 32 partial sums
    int lane = threadIdx.x % warpSize;
    int wid = threadIdx.x / warpSize;


    for (int offset = warpSize/2; offset > 0; offset /= 2) 
    {
        xMin = min(xMin, __shfl_xor_sync(FULL_MASK, xMin, offset));
        xMax = max(xMax, __shfl_xor_sync(FULL_MASK, xMax, offset));
        yMin = min(yMin, __shfl_xor_sync(FULL_MASK, yMin, offset));
        yMax = max(yMax, __shfl_xor_sync(FULL_MASK, yMax, offset));
        xSum += __shfl_down_sync(FULL_MASK, xSum, offset);
        x2Sum += __shfl_down_sync(FULL_MASK, x2Sum, offset);
        ySum += __shfl_down_sync(FULL_MASK, ySum, offset);
        y2Sum += __shfl_down_sync(FULL_MASK, y2Sum, offset);
        xySum += __shfl_down_sync(FULL_MASK, xySum, offset);
    }

    if (lane==0){
        shared[wid] = xMin;
        shared[32+wid] = xMax;
        shared[32*2+wid] = yMin;
        shared[32*3+wid] = yMax;
        shared[32*4+wid] = xSum;
        shared[32*5+wid] = x2Sum;
        shared[32*6+wid] = ySum;
        shared[32*7+wid] = y2Sum;
        shared[32*8+wid] = xySum;
        //printf("shared%i=%e\n", 32*6+wid, shared[32*6+wid]);
    }

    __syncthreads();                  

    if (threadIdx.x < blockDim.x / warpSize){
        xMin = shared[lane];
        xMax = shared[32+lane];
        yMin = shared[32*2+lane];
        yMax = shared[32*3+lane];
        xSum = shared[32*4+lane];
        x2Sum = shared[32*5+lane];
        ySum = shared[32*6+lane];
        y2Sum = shared[32*7+lane];
        xySum = shared[32*8+lane];
    }else{
        xMin = 0;  
        xMax = 0; 
        yMin = 0; 
        yMax = 0; 
        xSum = 0; 
        x2Sum = 0;
        ySum = 0; 
        y2Sum = 0;
        xySum = 0;
    }
    //if (wid==0)printf("ddata%i=%e:%e\n", 32*6+lane, shared[32*6+lane], ySum);

    if (wid==0){
        for (int offset = warpSize/2; offset > 0; offset /= 2) 
        {
            xMin = min(xMin, __shfl_xor_sync(FULL_MASK, xMin, offset));
            xMax = max(xMax, __shfl_xor_sync(FULL_MASK, xMax, offset));
            yMin = min(yMin, __shfl_xor_sync(FULL_MASK, yMin, offset));
            yMax = max(yMax, __shfl_xor_sync(FULL_MASK, yMax, offset));
            xSum += __shfl_down_sync(FULL_MASK, xSum, offset);
            x2Sum += __shfl_down_sync(FULL_MASK, x2Sum, offset);
            ySum += __shfl_down_sync(FULL_MASK, ySum, offset);
            y2Sum += __shfl_down_sync(FULL_MASK, y2Sum, offset);
            xySum += __shfl_down_sync(FULL_MASK, xySum, offset);
        }
    }
    double xMean=xSum/np;
    double yMean=ySum/np;
    double xSigma=sqrt(fabs((x2Sum/np)-(xMean*xMean)));
    double ySigma=sqrt(fabs((y2Sum/np)-(yMean*yMean)));
    double xyCov=(xySum/np)-(xMean*yMean);
    //if (wid==0)printf("ddata%i=%e:%e:%e\n", np, (xySum/np),(xMean*yMean), xyCov);

    double c1,c2;
    if(xMax-xMin==0){
      c1=K1*K1;
      c2=K2*K2;
    }else{
      c1=K1*K1*(xMax-xMin)*(xMax-xMin);
      c2=K2*K2*(xMax-xMin)*(xMax-xMin);
    }
    double c3=c2/2;
      
    double luminance=(2*xMean*yMean+c1)/(xMean*xMean+yMean*yMean+c1);
    double contrast=(2*xSigma*ySigma+c2)/(xSigma*xSigma+ySigma*ySigma+c2);
    double structure=(xyCov+c3)/(xSigma*ySigma+c3);
    double ssim=luminance*contrast*structure;
    //if (wid==0)printf("ddata%i=%e:%e:%e:%e:%e:%e\n", np,xMean,yMean,xSigma,ySigma,xyCov,ssim);
    return ssim;
    return ySum;
}

__global__ void ssim(float *data1, float *data2, double *results, int r3, int r2, int r1, int ssimSize, int ssimShift) 
{
    int tid = threadIdx.x;
    int bid = blockIdx.x;
    float xMin;  
    float xMax; 
    float yMin; 
    float yMax; 
    float xSum; 
    float x2Sum;
    float ySum; 
    float y2Sum;
    float xySum;
    int i;
    int np = ssimSize * ssimSize * ssimSize;
    double result = 0;

    int Offsetx;
    int Offsety;
    int Offsetz = bid * ssimShift;
    //int Offsetz = 0;
    int x=0;

    for (Offsetx=0; Offsetx+ssimSize<=r3; Offsetx+=ssimShift){
        for (Offsety=0; Offsety+ssimSize<=r2; Offsety+=ssimShift){
    //for (Offsetx=0; Offsetx+ssimSize<=ssimSize; Offsetx+=ssimShift){
    //    for (Offsety=0; Offsety+ssimSize<=ssimSize; Offsety+=ssimShift){
    //for (Offsetx=36; Offsetx<=36; Offsetx+=ssimShift){
    //    for (Offsety=330; Offsety<=331; Offsety+=ssimShift){
            xMin = data1[0];  
            xMax = data1[0]; 
            yMin = data2[0]; 
            yMax = data2[0]; 
            xSum = 0; 
            x2Sum = 0;
            ySum = 0; 
            y2Sum = 0;
            xySum = 0;

            for (i=tid; i<ssimSize*ssimSize*ssimSize; i+=blockDim.x)
            {
                int Winx = i / (ssimSize*ssimSize);
                int Winy = i % (ssimSize*ssimSize) / ssimSize;
                int Winz = i % (ssimSize*ssimSize) % ssimSize;
                int index = (Offsetx + Winx) * r1 * r2 + (Offsety + Winy) * r2 + (Offsetz + Winz);
                float xdata = data1[index];
                float ydata = data2[index];
                //if (i<blockDim.x)
                //{
                //    xMin = xdata;  
                //    xMax = xdata; 
                //    yMin = ydata; 
                //    yMax = ydata; 
                //}
                xMin = min(xMin, xdata);
                xMax = max(xMax, xdata);
                yMin = min(yMin, ydata);
                yMax = max(yMax, ydata);
                xSum += xdata;
                x2Sum += xdata * xdata;
                ySum += ydata;
                y2Sum += ydata * ydata;
                xySum += xdata * ydata;
            //printf("ydata%i=%e\n",index, xdata);

            }
            //printf("ydata=%e\n", ySum);
            //if (tid == 0) results[bid] = blockReduceSum(val);
            result += ssim_winComp(xMin, xMax, yMin, yMax, xSum, x2Sum, ySum, y2Sum, xySum, np);
            //results[x++] = ssim_winComp(xMin, xMax, yMin, yMax, xSum, x2Sum, ySum, y2Sum, xySum, np);
            //results[tid] = blockReduceSum(val);
            __syncthreads();                  

        }
    }
if (tid == 0)printf("ydata%i,%i=%e\n",Offsetx,Offsety, result);
    if (tid==0) results[bid] = result;
    
}
