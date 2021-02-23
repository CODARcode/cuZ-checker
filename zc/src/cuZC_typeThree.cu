#include <stdio.h>
#include <math.h>
#include "cuZC_ssim.h"
#include "cuZC_typeThree.h"
#include "matrix.hpp"
#include <cooperative_groups.h>

namespace cg = cooperative_groups;

__global__ void type_three(float *data1, float *data2, double *results, int r3, int r2, int r1, int ssimSize, int ssimShift, int yNum) 
{
    //cg::grid_group grid = cg::this_grid();
    float xMin;  
    float xMax; 
    float yMin; 
    float yMax; 
    float xSum; 
    float x2Sum;
    float ySum; 
    float y2Sum;
    float xySum;
    float xdata, ydata, xdata_shflx, ydata_shflx, xdata_shfld, ydata_shfld;

    int i, j;
    int xNum = (warpSize-ssimSize)/ssimShift+1;
    int h=blockIdx.x*yNum*ssimShift;
    if ((h+blockDim.y)>r2) yNum = (r2-h-ssimSize)/ssimShift+1;
    int wstride=xNum*ssimShift;
    double np = ssimSize * ssimSize * ssimSize;
    unsigned mask;
    
    int wsize = (r1-ssimSize+ssimShift)/wstride + ((r1-ssimSize+ssimShift)%wstride?1:0);
    
    //static __shared__ float shared[9*xNum*(yNum*ssimSize+blockDim.y)];
    //static __shared__ float shared[9*26*(2*7+8)];
    extern __shared__ float shared[];

    for (int w=0; w<wsize*wstride; w+=wstride){
        if ((w+blockDim.x)>r1) xNum = (r1-w-ssimSize)/ssimShift+1;
        if (w+threadIdx.x<r1){
            for (int l=0; l<r3; l++){
            //for (int l=0; l<ssimSize; l++){

                xMin = data1[0];  
                xMax = data1[0]; 
                yMin = data2[0]; 
                yMax = data2[0]; 
                xdata= data1[l*r1*r2+(h+threadIdx.y)*r1+(w+threadIdx.x)];
                ydata= data2[l*r1*r2+(h+threadIdx.y)*r1+(w+threadIdx.x)];
                xMin = min(xMin, xdata);
                xMax = max(xMax, xdata);
                yMin = min(yMin, ydata);
                yMax = max(yMax, ydata);
                xSum = xdata; 
                x2Sum = xdata * xdata;
                ySum = ydata; 
                y2Sum = ydata * ydata;
                xySum = xdata * ydata;

                mask = __ballot_sync(FULL_MASK, threadIdx.x < (r1-w));

                for (int offset = 1; offset<ssimSize; offset++) 
                {
                    //xdata_shflx = __shfl_xor_sync(mask, xdata, offset);
                    //ydata_shflx = __shfl_xor_sync(mask, ydata, offset);
                    xdata_shfld = __shfl_down_sync(mask, xdata, offset);
                    ydata_shfld = __shfl_down_sync(mask, ydata, offset);

                    xMin = min(xMin, xdata_shfld);
                    xMax = max(xMax, xdata_shfld);
                    yMin = min(yMin, ydata_shfld);
                    yMax = max(yMax, ydata_shfld);
                    xSum += xdata_shfld;
                    x2Sum += xdata_shfld * xdata_shfld;
                    ySum += ydata_shfld;
                    y2Sum += ydata_shfld * ydata_shfld;
                    xySum += xdata_shfld * ydata_shfld;
                }
                
                if (threadIdx.x%ssimShift==0 && threadIdx.x/ssimShift<xNum){
//if (xMax!=0.0) printf("test%i,%i,%i=%e\n", h,w,l,xMax);
                    shared[(9*yNum*ssimSize+threadIdx.y)*xNum+threadIdx.x/ssimShift] = xMin;
                    shared[(9*yNum*ssimSize+blockDim.y+threadIdx.y)*xNum+threadIdx.x/ssimShift] = yMin;
                    shared[(9*yNum*ssimSize+2*blockDim.y+threadIdx.y)*xNum+threadIdx.x/ssimShift] = xMax;
                    shared[(9*yNum*ssimSize+3*blockDim.y+threadIdx.y)*xNum+threadIdx.x/ssimShift] = yMax;
                    shared[(9*yNum*ssimSize+4*blockDim.y+threadIdx.y)*xNum+threadIdx.x/ssimShift] = xSum;
                    shared[(9*yNum*ssimSize+5*blockDim.y+threadIdx.y)*xNum+threadIdx.x/ssimShift] = x2Sum;
                    shared[(9*yNum*ssimSize+6*blockDim.y+threadIdx.y)*xNum+threadIdx.x/ssimShift] = ySum;
                    shared[(9*yNum*ssimSize+7*blockDim.y+threadIdx.y)*xNum+threadIdx.x/ssimShift] = y2Sum;
                    shared[(9*yNum*ssimSize+8*blockDim.y+threadIdx.y)*xNum+threadIdx.x/ssimShift] = xySum;

                }
                //if (threadIdx.x==0) printf("test%i,%i=%e\n",l,threadIdx.y,ySum);
                __syncthreads();                  

                if (threadIdx.x<xNum){
                    for (j=0;j<yNum;j++){
                        if (threadIdx.y==j){
                            for (i=j;i<(ssimSize+j);i++) xMin = min(xMin, shared[(9*yNum*ssimSize+i)*xNum+threadIdx.x]);
                            shared[9*(yNum*(l%ssimSize)+j)*xNum+threadIdx.x] = xMin;
                        }else if (threadIdx.y==(j+1)){
                            for (i=j;i<(ssimSize+j);i++) yMin = min(yMin, shared[(9*yNum*ssimSize+blockDim.y+i)*xNum+threadIdx.x]);
                            shared[(9*(yNum*(l%ssimSize)+j)+1)*xNum+threadIdx.x] = yMin;
                        }else if (threadIdx.y==(j+2)){
                            //xMax = shared[(9*yNum*ssimSize+2*blockDim.y)*xNum+threadIdx.x];
                            for (i=j;i<(ssimSize+j);i++) xMax = max(xMax, shared[(9*yNum*ssimSize+2*blockDim.y+i)*xNum+threadIdx.x]);
                            shared[(9*(yNum*(l%ssimSize)+j)+2)*xNum+threadIdx.x] = xMax;
                        }else if (threadIdx.y==(j+3)){
                            for (i=j;i<(ssimSize+j);i++) yMax = max(yMax, shared[(9*yNum*ssimSize+3*blockDim.y+i)*xNum+threadIdx.x]);
                            shared[(9*(yNum*(l%ssimSize)+j)+3)*xNum+threadIdx.x] = yMax;
                        }

                        xSum = 0;
                        if (threadIdx.y<5){
                            for (i=j;i<(ssimSize+j);i++) xSum += shared[(9*yNum*ssimSize+(4+threadIdx.y)*blockDim.y+i)*xNum+threadIdx.x];
                            shared[(9*(yNum*(l%ssimSize)+j)+(4+threadIdx.y))*xNum+threadIdx.x] = xSum;
                        }
                    }
                }
                __syncthreads();                  

                if (l>(ssimSize-2)){
                    if ((l-ssimSize+1)%ssimShift==0){
                        if (threadIdx.x<xNum){
                            if (threadIdx.y<yNum){
                                xMin = shared[(9*(yNum*0+threadIdx.y)+0)*xNum+threadIdx.x];  
                                yMin = shared[(9*(yNum*0+threadIdx.y)+1)*xNum+threadIdx.x]; 
                                xMax = shared[(9*(yNum*0+threadIdx.y)+2)*xNum+threadIdx.x]; 
                                yMax = shared[(9*(yNum*0+threadIdx.y)+3)*xNum+threadIdx.x]; 
                                xSum = shared[(9*(yNum*0+threadIdx.y)+4)*xNum+threadIdx.x]; 
                                x2Sum =shared[(9*(yNum*0+threadIdx.y)+5)*xNum+threadIdx.x];
                                ySum = shared[(9*(yNum*0+threadIdx.y)+6)*xNum+threadIdx.x]; 
                                y2Sum =shared[(9*(yNum*0+threadIdx.y)+7)*xNum+threadIdx.x];
                                xySum =shared[(9*(yNum*0+threadIdx.y)+8)*xNum+threadIdx.x];

                                for (i=1;i<ssimSize;i++) {
                                    xMin = min(xMin, shared[(9*(yNum*i+threadIdx.y)+0)*xNum+threadIdx.x]);
                                    yMin = min(yMin, shared[(9*(yNum*i+threadIdx.y)+1)*xNum+threadIdx.x]);
                                    xMax = max(xMax, shared[(9*(yNum*i+threadIdx.y)+2)*xNum+threadIdx.x]);
                                    yMax = max(yMax, shared[(9*(yNum*i+threadIdx.y)+3)*xNum+threadIdx.x]);
                                    xSum += shared[(9*(yNum*i+threadIdx.y)+4)*xNum+threadIdx.x];
                                    x2Sum += shared[(9*(yNum*i+threadIdx.y)+5)*xNum+threadIdx.x];
                                    ySum += shared[(9*(yNum*i+threadIdx.y)+6)*xNum+threadIdx.x];
                                    y2Sum += shared[(9*(yNum*i+threadIdx.y)+7)*xNum+threadIdx.x];
                                    xySum += shared[(9*(yNum*i+threadIdx.y)+8)*xNum+threadIdx.x];

                                }
                //if (threadIdx.x==0) printf("test%i,%i=%e\n",l,(l-ssimSize+1)/ssimShift*30*yNum+30*threadIdx.y+threadIdx.x,ySum);

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

                                results[(l-ssimSize+1)/ssimShift*r1*r2+r1*(h+threadIdx.y)+(w+threadIdx.x)] = ssim;
                            }
                        }
                        __syncthreads();                  
                    }
                }
            }

        }

    }
    //cg::sync(grid);
}

__global__ void gridR_typeThree(double *results, int size) 
{
    int tidx = threadIdx.x;
    double data = results[tidx];

    for (int i=(tidx+blockDim.x); i<size; i+=blockDim.x)
        data += results[i];
    __syncthreads();                  

    for (int offset = warpSize/2; offset > 0; offset /= 2) 
        data += __shfl_down_sync(FULL_MASK, data, offset);

    if (tidx==0) results[0] = data;
        
}
