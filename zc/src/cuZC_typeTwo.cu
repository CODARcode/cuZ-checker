#include <stdio.h>
#include <math.h>
#include "cuZC_ssim.h"
#include "cuZC_typeTwo.h"
#include "matrix.hpp"

__global__ void type_two(float *data, float *der, float *autocor, int r3, int r2, int r1, float avg, size_t order) 
{
    int tidx = threadIdx.x;
    int tidy = threadIdx.y;

    float base, sum;
    int i, j;
    int h=blockIdx.x*(16-order*2);
    double dx, dy, dz;
    
    int wsize = (r2-order*2)/(16-order*2) + ((r2-order*2)%(16-order*2)?1:0);
    int lsize = (r1-order*2)/(16-order*2) + ((r1-order*2)%(16-order*2)?1:0);
    
    //static __shared__ float shared[16*16*16];
    extern __shared__ float shared[];
    float *bdata = shared;
    float *cor = &shared[blockDim.x * blockDim.y * 16];
    unsigned mask;

    for (int w=0; w<wsize*(16-order*2); w+=(16-order*2)){
        for (int l=0; l<lsize*(16-order*2); l+=(16-order*2)){
            for (i=0; i<16; i++){
                if ((h+i)<r3 && (l+tidx)<r1 && (w+tidy)<r2){
                    bdata[i*16*16+tidy*16+tidx] = data[(h+i)*r1*r2+(w+tidy)*r1+(l+tidx)];
                }
                    
            }
            __syncthreads();                  

            for (i=0; i<(16-order*2); i++){
                if (tidx<(16-order*2) && tidy<(16-order*2)){
                    if ((h+i)<(r3-order*2) && (l+tidx)<(r1-order*2) && (w+tidy)<(r2-order*2)){
                        base = bdata[(i+order)*16*16+(tidy+order)*16+tidx];
                        dx = (bdata[(i+order)*16*16+(tidy+order)*16+tidx+order*2] - base)/2;
                        base = bdata[(i+order)*16*16+tidy*16+tidx+order];
                        dy = (bdata[(i+order)*16*16+(tidy+order*2)*16+tidx+order] - base)/2;
                        base = bdata[i*16*16+(tidy+order)*16+tidx+order];
                        dz = (bdata[(i+order*2)*16*16+(tidy+order)*16+tidx+order] - base)/2;
                        //if (bid==0)printf("index=%i,%e,%i,%i\n",i+2,dz,h+tidy,l+tidx);
                        //if (Data!=0.0) printf("ddata%i,%i,%i,%i,%i,%i=%e\n",w,l,bid,i,tidx,tidy,sqrt(dx*dx+dy*dy+dz*dz));
                        der[(h+i)*(r1-order*2)*(r2-order*2)+(w+tidy)*(r1-order*2)+(l+tidx)] = sqrt(dx*dx+dy*dy+dz*dz);
                        //if (der[(w+i)*(r1-order*2)*(r2-order*2)+(h+tidy)*(r1-order*2)+(l+tidx)]!=0.0) printf("ddata%i=%e\n",(w+i)*(r1-order*2)*(r2-order*2)+(h+tidy)*(r1-order*2)+(l+tidx),der[(w+i)*(r1-order*2)*(r2-order*2)+(h+tidy)*(r1-order*2)+(l+tidx)]);

                        mask = __ballot_sync(FULL_MASK, 1);
                        base = bdata[i*16*16+tidy*16+tidx];

                        for (j=1; j<=order*2; j++){
                            sum = (bdata[(i+j)*16*16+(tidy+j)*16+tidx+j]-avg) * (base-avg);

                            for (int offset = warpSize/2; offset > 0; offset /= 2) 
                                sum += __shfl_down_sync(mask, sum, offset);

                            if (tidx==0) cor[blockDim.y*(j-1)+tidy] = sum;
                        }
                    }
                }
                __syncthreads();                  

                if (tidy<order*2){
                    if (tidx < (16-order*2) && (w+tidx)<(r2-order*2))
                    {
                        sum = cor[blockDim.y*tidy+tidx];
                        mask = __ballot_sync(FULL_MASK, 1);
                    } else sum = 0;
                    for (int offset = warpSize/2; offset > 0; offset /= 2) 
                        sum += __shfl_down_sync(mask, sum, offset);

                    if (tidx==0) autocor[gridDim.x*tidy+blockIdx.x] += sum;
                }
                __syncthreads();                  
            }
        }
    }
}
