#include <stdio.h>
#include <math.h>
#include "cuZC_ssim.h"
#include "cuZC_typeThree.h"
#include "matrix.hpp"

__global__ void type_three(float *data1, float *data2, double *results, int r3, int r2, int r1, int ssimSize, int ssimShift) 
{
    int tidx = threadIdx.x;
    int tidy = threadIdx.y;
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
    float xdata, ydata, xdata_shflx, ydata_shflx, xdata_shfld, ydata_shfld;

    int i, j;
    int xNum = (warpSize-ssimSize)/ssimShift+1;
    int yNum = (blockDim.y-ssimSize)/ssimShift+1;
    int h=bid*yNum*ssimShift;
    int wstride=xNum*ssimShift;
    double np = ssimSize * ssimSize * ssimSize;
    
    int wsize = (r1-ssimSize+ssimShift)/wstride + ((r1-ssimSize+ssimShift)%wstride?1:0);
    
    //static __shared__ float shared[9*xNum*(yNum*ssimSize+blockDim.y)];
    static __shared__ float shared[9*26*(2*7+8)];

    //for (int w=0; w<wsize*wstride; w+=wstride){
    int w=0;
        for (int l=0; l<r3; l++){
            xMin = data1[0];  
            xMax = data1[0]; 
            yMin = data2[0]; 
            yMax = data2[0]; 
            xdata= data1[l*r1*r2+(h+tidy)*r1+(w+tidx)];
            ydata= data2[l*r1*r2+(h+tidy)*r1+(w+tidx)];
            xMin = min(xMin, xdata);
            xMax = max(xMax, xdata);
            yMin = min(yMin, ydata);
            yMax = max(yMax, ydata);
            xSum = xdata; 
            x2Sum = xdata * xdata;
            ySum = ydata; 
            y2Sum = ydata * ydata;
            xySum = xdata * ydata;


            for (int offset = 1; offset<ssimSize; offset++) 
            {
                xdata_shflx = __shfl_xor_sync(FULL_MASK, xdata, offset);
                ydata_shflx = __shfl_xor_sync(FULL_MASK, ydata, offset);
                xdata_shfld = __shfl_down_sync(FULL_MASK, xdata, offset);
                ydata_shfld = __shfl_down_sync(FULL_MASK, ydata, offset);

                xMin = min(xMin, xdata_shflx);
                xMax = max(xMax, xdata_shflx);
                yMin = min(yMin, ydata_shflx);
                yMax = max(yMax, ydata_shflx);
                xSum += xdata_shfld;
                x2Sum += xdata_shfld * xdata_shfld;
                ySum += ydata_shfld;
                y2Sum += ydata_shfld * ydata_shfld;
                xySum += xdata_shfld * ydata_shfld;
            }
            
            if (tidx<xNum){
                shared[(9*yNum*ssimSize+tidy)*xNum+tidx] = xMin;
                shared[(9*yNum*ssimSize+blockDim.y+tidy)*xNum+tidx] = yMin;
                shared[(9*yNum*ssimSize+2*blockDim.y+tidy)*xNum+tidx] = xMax;
                shared[(9*yNum*ssimSize+3*blockDim.y+tidy)*xNum+tidx] = yMax;
                shared[(9*yNum*ssimSize+4*blockDim.y+tidy)*xNum+tidx] = xSum;
                shared[(9*yNum*ssimSize+5*blockDim.y+tidy)*xNum+tidx] = x2Sum;
                shared[(9*yNum*ssimSize+6*blockDim.y+tidy)*xNum+tidx] = ySum;
                shared[(9*yNum*ssimSize+7*blockDim.y+tidy)*xNum+tidx] = y2Sum;
                shared[(9*yNum*ssimSize+8*blockDim.y+tidy)*xNum+tidx] = xySum;

            }
            //if (tidx==0) printf("test%i,%i=%e\n",l,tidy,ySum);
            __syncthreads();                  

            if (tidx<xNum){
                for (j=0;j<yNum;j++){
                    if (tidy==j){
                        for (i=j;i<(ssimSize+j);i++) xMin = min(xMin, shared[(9*yNum*ssimSize+i)*xNum+tidx]);
                        shared[9*(yNum*(l%ssimSize)+j)*xNum+tidx] = xMin;
                    }else if (tidy==(j+1)){
                        for (i=j;i<(ssimSize+j);i++) yMin = min(yMin, shared[(9*yNum*ssimSize+blockDim.y+i)*xNum+tidx]);
                        shared[(9*(yNum*(l%ssimSize)+j)+1)*xNum+tidx] = yMin;
                    }else if (tidy==(j+2)){
                        for (i=j;i<(ssimSize+j);i++) xMax = max(xMax, shared[(9*yNum*ssimSize+2*blockDim.y+i)*xNum+tidx]);
                        shared[(9*(yNum*(l%ssimSize)+j)+2)*xNum+tidx] = xMax;
                    }else if (tidy==(j+3)){
                        for (i=j;i<(ssimSize+j);i++) yMax = max(yMax, shared[(9*yNum*ssimSize+3*blockDim.y+i)*xNum+tidx]);
                        shared[(9*(yNum*(l%ssimSize)+j)+3)*xNum+tidx] = yMax;
                    }

                    xSum = 0;
                    if (tidy<5){
                        for (i=j;i<(ssimSize+j);i++) xSum += shared[(9*yNum*ssimSize+(4+tidy)*blockDim.y+i)*xNum+tidx];
                        shared[(9*(yNum*(l%ssimSize)+j)+(4+tidy))*xNum+tidx] = xSum;
                    }
                }
            }
            __syncthreads();                  

            if (l>(ssimSize-2)){
                if ((l-ssimSize+1)%ssimShift==0){
                    if (tidx<xNum){
                        if (tidy<yNum){
                            xMin = shared[(9*(yNum*0+tidy)+0)*xNum+tidx];  
                            yMin = shared[(9*(yNum*0+tidy)+1)*xNum+tidx]; 
                            xMax = shared[(9*(yNum*0+tidy)+2)*xNum+tidx]; 
                            yMax = shared[(9*(yNum*0+tidy)+3)*xNum+tidx]; 
                            xSum = shared[(9*(yNum*0+tidy)+4)*xNum+tidx]; 
                            x2Sum =shared[(9*(yNum*0+tidy)+5)*xNum+tidx];
                            ySum = shared[(9*(yNum*0+tidy)+6)*xNum+tidx]; 
                            y2Sum =shared[(9*(yNum*0+tidy)+7)*xNum+tidx];
                            xySum =shared[(9*(yNum*0+tidy)+8)*xNum+tidx];

                            for (i=1;i<ssimSize;i++) {
                                xMin = min(xMin, shared[(9*(yNum*i+tidy)+0)*xNum+tidx]);
                                yMin = min(yMin, shared[(9*(yNum*i+tidy)+1)*xNum+tidx]);
                                xMax = max(xMax, shared[(9*(yNum*i+tidy)+2)*xNum+tidx]);
                                yMax = max(yMax, shared[(9*(yNum*i+tidy)+3)*xNum+tidx]);
                                xSum += shared[(9*(yNum*i+tidy)+4)*xNum+tidx];
                                x2Sum += shared[(9*(yNum*i+tidy)+5)*xNum+tidx];
                                ySum += shared[(9*(yNum*i+tidy)+6)*xNum+tidx];
                                y2Sum += shared[(9*(yNum*i+tidy)+7)*xNum+tidx];
                                xySum += shared[(9*(yNum*i+tidy)+8)*xNum+tidx];

                            }
            if (tidx==0) printf("test%i,%i=%e\n",l,(l-ssimSize+1)/ssimShift*30*yNum+30*tidy+tidx,ySum);

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

                            results[(l-ssimSize+1)/ssimShift*30*yNum+30*tidy+tidx] = ySum;
                        }
                    }
                    __syncthreads();                  
                }
            }
        }

    //}
}
