#include <stdio.h>
#include <math.h>
#include "cuZC_ssim.h"
#include "cuZC_typeTwo.h"
#include "matrix.hpp"

__global__ void type_two(float *data, float *der, int r3, int r2, int r1, size_t order) 
{
    int tidx = threadIdx.x;
    int tidy = threadIdx.y;
    int bid = blockIdx.x;

    float Data;
    int i, j;
    int h=bid*(16-order*2);
    double dx, dy, dz;
    
    int wsize = (r2-order*2)/(16-order*2) + ((r2-order*2)%(16-order*2)?1:0);
    int lsize = (r1-order*2)/(16-order*2) + ((r1-order*2)%(16-order*2)?1:0);
    
    static __shared__ float shared[16*16*16];

    for (int w=0; w<wsize*(16-order*2); w+=(16-order*2)){
        for (int l=0; l<lsize*(16-order*2); l+=(16-order*2)){
            for (i=0; i<16; i++){
                if ((h+i)<r3 && (l+tidx)<r1 && (w+tidy)<r2){
                    shared[i*16*16+tidy*16+tidx] = data[(h+i)*r1*r2+(w+tidy)*r1+(l+tidx)];
                }
                    
            }
            __syncthreads();                  

            for (i=0; i<(16-order*2); i++){
                if (tidx<(16-order*2) && tidy<(16-order*2)){
                    if ((h+i)<(r3-order*2) && (l+tidx)<(r1-order*2) && (w+tidy)<(r2-order*2)){
                        Data = shared[i*16*16+tidy*16+tidx];
                        dx = (shared[i*16*16+tidy*16+tidx+order*2] - Data)/2;
                        dy = (shared[i*16*16+(tidy+order*2)*16+tidx] - Data)/2;
                        dz = (shared[(i+order*2)*16*16+tidy*16+tidx] - Data)/2;
                        //if (bid==0)printf("index=%i,%e,%i,%i\n",i+2,dz,h+tidy,l+tidx);
                        //if (Data!=0.0) printf("ddata%i,%i,%i,%i,%i,%i=%e\n",w,l,bid,i,tidx,tidy,sqrt(dx*dx+dy*dy+dz*dz));
                        der[(h+i)*(r1-order*2)*(r2-order*2)+(w+tidy)*(r1-order*2)+(l+tidx)] = sqrt(dx*dx+dy*dy+dz*dz);
                        //der[(h+i)*(r1-order*2)*(r2-order*2)+(w+tidy)*(r1-order*2)+(l+tidx)] = Data;
                        //if (der[(w+i)*(r1-order*2)*(r2-order*2)+(h+tidy)*(r1-order*2)+(l+tidx)]!=0.0) printf("ddata%i=%e\n",(w+i)*(r1-order*2)*(r2-order*2)+(h+tidy)*(r1-order*2)+(l+tidx),der[(w+i)*(r1-order*2)*(r2-order*2)+(h+tidy)*(r1-order*2)+(l+tidx)]);
                    
                    }
                }
            }
            __syncthreads();                  
        }
    }

}

