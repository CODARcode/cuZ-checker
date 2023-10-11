#include<cmath>
#include"cuZC_derivatives.h"

// __global__ void testKernel(float* data,int nx,int ny,int nz){
//     int n=nx*ny*nz;
//     float maxval=-inf;
//     for(int i=0;i<n;i++) maxval=max(maxval,data[i]);
//     printf("+++%.12f\n",maxval);
// }
//#ifndef float
//#define float float
//#endif

#ifndef arr3
#define arr3(x, y, z) ((x)*ny*nz+(y)*nz+(z))
#endif

#ifndef square
#define square(x) ((x)*(x))
#endif

#ifndef _BLOCK_SZ
#define BLOCKSZX 8
#define BLOCKSZY 8
#define BLOCKSZZ 8
#endif

#define blksz 1024
// must be 2^k, and 0<= k <= 10

#ifndef _INF
#define _INF
const float inf = 1. / 0.;
#endif


__global__ void maximumReduction(float *data, int p, int n) {

    int id = (blockIdx.x * blockDim.x + threadIdx.x) * p;

    if (id >= n) return;

    for (int stride = p; stride < p * blksz; stride <<= 1) {

        if (id % (stride << 1) == 0) {

            if (id + stride < n) {

                if (data[id + stride] > data[id]) {

                    data[id] = data[id + stride];
                }
            }
        }

        __syncthreads();
    }
}

float findMaximumUsingReduction(float *data, int n) {

    float *tem;

    cudaMalloc(&tem, n * sizeof(float));

    cudaMemcpy(tem, data, n * sizeof(float), cudaMemcpyDeviceToDevice);

    int p = 1;

    while (true) {

        int blknum = n / (p * blksz) + (n % (p * blksz) > 0);

        maximumReduction<<<blknum, blksz>>>(tem, p, n);

        cudaDeviceSynchronize();

        if (blknum == 1) break;

        p *= blksz;
    }

    float ans;

    cudaMemcpy(&ans, tem, sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(tem);

    return ans;
}

__global__ void minimumReduction(float *data, int p, int n) {

    int id = (blockIdx.x * blockDim.x + threadIdx.x) * p;

    if (id >= n) return;

    for (int stride = p; stride < p * blksz; stride <<= 1) {

        if (id % (stride << 1) == 0) {

            if (id + stride < n) {

                if (data[id + stride] < data[id]) {

                    data[id] = data[id + stride];
                }
            }
        }

        __syncthreads();
    }
}

float findMinimumUsingReduction(float *data, int n) {

    float *tem;

    cudaMalloc(&tem, n * sizeof(float));

    cudaMemcpy(tem, data, n * sizeof(float), cudaMemcpyDeviceToDevice);

    int p = 1;

    while (true) {

        int blknum = n / (p * blksz) + (n % (p * blksz) > 0);

        minimumReduction<<<blknum, blksz>>>(tem, p, n);

        cudaDeviceSynchronize();

        if (blknum == 1) break;

        p *= blksz;
    }

    float ans;

    cudaMemcpy(&ans, tem, sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(tem);

    return ans;
}

__global__ void sumReduction(float *data, int p, int n) {

    int id = (blockIdx.x * blockDim.x + threadIdx.x) * p;

    if (id >= n) return;

    for (int stride = p; stride < p * blksz; stride <<= 1) {

        if (id % (stride << 1) == 0) {

            if (id + stride < n) {

                data[id] += data[id + stride];
            }
        }

        __syncthreads();
    }
}

float sumupUsingReduction(float *data, int n) {

    float *tem;

    cudaMalloc(&tem, n * sizeof(float));

    cudaMemcpy(tem, data, n * sizeof(float), cudaMemcpyDeviceToDevice);

    int p = 1;

    while (true) {

        int blknum = n / (p * blksz) + (n % (p * blksz) > 0);

        sumReduction<<<blknum, blksz>>>(tem, p, n);

        cudaDeviceSynchronize();

        if (blknum == 1) break;

        p *= blksz;
    }

    float ans;

    cudaMemcpy(&ans, tem, sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(tem);

    return ans;
}

__global__ void sumReductionForceDouble(double *data, int p, int n) {

    int id = (blockIdx.x * blockDim.x + threadIdx.x) * p;

    if (id >= n) return;

    for (int stride = p; stride < p * blksz; stride <<= 1) {

        if (id % (stride << 1) == 0) {

            if (id + stride < n) {

                data[id] += data[id + stride];
            }
        }

        __syncthreads();
    }
}

double sumupUsingReductionForceDouble(double *data, int n) {

    double *tem;

    cudaMalloc(&tem, n * sizeof(double));

    cudaMemcpy(tem, data, n * sizeof(double), cudaMemcpyDeviceToDevice);

    int p = 1;

    while (true) {

        int blknum = n / (p * blksz) + (n % (p * blksz) > 0);

        sumReductionForceDouble<<<blknum, blksz>>>(tem, p, n);

        cudaDeviceSynchronize();

        if (blknum == 1) break;

        p *= blksz;
    }

    double ans;

    cudaMemcpy(&ans, tem, sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(tem);

    return ans;
}

__global__ void
derivativesKernel(float *data, float *dx, float *dy, float *dz, float *dx2, float *dy2, float *dz2, float *dxy, float *dyz, float *dzx, float *gl,
                  float *lap, int nx, int ny, int nz) {

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x < nx && y < ny && z < nz) {

        // 1st order derivatives

        if (x < nx - 1) {
            dx[arr3(x, y, z)] = data[arr3(x + 1, y, z)] - data[arr3(x, y, z)];
        } else {
            dx[arr3(x, y, z)] = 0;
        }

        if (y < ny - 1) {
            dy[arr3(x, y, z)] = data[arr3(x, y + 1, z)] - data[arr3(x, y, z)];
        } else {
            dy[arr3(x, y, z)] = 0;
        }

        if (z < nz - 1) {
            dz[arr3(x, y, z)] = data[arr3(x, y, z + 1)] - data[arr3(x, y, z)];
        } else {
            dz[arr3(x, y, z)] = 0;
        }

        // 2rd order derivatives

        if (x < nx - 2) {
            dx2[arr3(x, y, z)] = data[arr3(x + 2, y, z)] - 2. * data[arr3(x + 1, y, z)] + data[arr3(x, y, z)];
        } else {
            dx2[arr3(x, y, z)] = 0;
        }

        if (y < ny - 2) {
            dy2[arr3(x, y, z)] = data[arr3(x, y + 2, z)] - 2. * data[arr3(x, y + 1, z)] + data[arr3(x, y, z)];
        } else {
            dy2[arr3(x, y, z)] = 0;
        }

        if (z < nz - 2) {
            dz2[arr3(x, y, z)] = data[arr3(x, y, z + 2)] - 2. * data[arr3(x, y, z + 1)] + data[arr3(x, y, z)];
        } else {
            dz2[arr3(x, y, z)] = 0;
        }

        if (x < nx - 1 && y < ny - 1) {

            dxy[arr3(x, y, z)] = data[arr3(x + 1, y + 1, z)] - data[arr3(x + 1, y, z)] - data[arr3(x, y + 1, z)] + data[arr3(x, y, z)];
        } else {

            dxy[arr3(x, y, z)] = 0;
        }

        if (y < ny - 1 && z < nz - 1) {

            dyz[arr3(x, y, z)] = data[arr3(x, y + 1, z + 1)] - data[arr3(x, y + 1, z)] - data[arr3(x, y, z + 1)] + data[arr3(x, y, z)];
        } else {

            dyz[arr3(x, y, z)] = 0;
        }

        if (z < nz - 1 && x < nx - 1) {

            dzx[arr3(x, y, z)] = data[arr3(x + 1, y, z + 1)] - data[arr3(x, y, z + 1)] - data[arr3(x + 1, y, z)] + data[arr3(x, y, z)];
        } else {

            dzx[arr3(x, y, z)] = 0;
        }

        //gradient length

        gl[arr3(x, y, z)] = sqrt(square(dx[arr3(x, y, z)]) + square(dy[arr3(x, y, z)]) + square(dz[arr3(x, y, z)]));

        //laplacian

        lap[arr3(x, y, z)] = square(dx2[arr3(x, y, z)]) + square(dy2[arr3(x, y, z)]) + square(dz2[arr3(x, y, z)]);
    }
}

__global__ void
squareMatrixKernelUsingBySobolev(float *data, float *dx, float *dy, float *dz, float *dx2, float *dy2, float *dz2, float *dxy, float *dyz, float *dzx,
                                 int nx, int ny, int nz) {

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x < nx && y < ny && z < nz) {

        data[arr3(x, y, z)] *= data[arr3(x, y, z)];
        dx[arr3(x, y, z)] *= dx[arr3(x, y, z)];
        dy[arr3(x, y, z)] *= dy[arr3(x, y, z)];
        dz[arr3(x, y, z)] *= dz[arr3(x, y, z)];
        dx2[arr3(x, y, z)] *= dx2[arr3(x, y, z)];
        dy2[arr3(x, y, z)] *= dy2[arr3(x, y, z)];
        dz2[arr3(x, y, z)] *= dz2[arr3(x, y, z)];
        dxy[arr3(x, y, z)] *= dxy[arr3(x, y, z)];
        dyz[arr3(x, y, z)] *= dyz[arr3(x, y, z)];
        dzx[arr3(x, y, z)] *= dzx[arr3(x, y, z)];

    }
}

std::vector<float>
sobolev(float *data, float *dx, float *dy, float *dz, float *dx2, float *dy2, float *dz2, float *dxy, float *dyz, float *dzx, int nx, int ny,
        int nz) {

    dim3 blocksz = dim3(BLOCKSZX, BLOCKSZY, BLOCKSZZ);
    dim3 blocknum = dim3(nx / BLOCKSZX + (nx % BLOCKSZX > 0), ny / BLOCKSZY + (ny % BLOCKSZY > 0), nz / BLOCKSZZ + (nz % BLOCKSZZ > 0));

    squareMatrixKernelUsingBySobolev<<<blocknum, blocksz>>>(data, dx, dy, dz, dx2, dy2, dz2, dxy, dyz, dzx, nx, ny, nz);

    cudaDeviceSynchronize();

    int n = nx * ny * nz;

    std::vector<float> s(3);

    s[0] = 0;
    s[0] += sumupUsingReduction(data, n);
    s[1] = s[0];
    s[1] += sumupUsingReduction(dx, n);
    s[1] += sumupUsingReduction(dy, n);
    s[1] += sumupUsingReduction(dz, n);
    s[2] = s[1];
    s[2] += sumupUsingReduction(dx2, n);
    s[2] += sumupUsingReduction(dy2, n);
    s[2] += sumupUsingReduction(dz2, n);
    s[2] += sumupUsingReduction(dxy, n);
    s[2] += sumupUsingReduction(dyz, n);
    s[2] += sumupUsingReduction(dzx, n);

    s[0] = sqrt(s[0] / n);
    s[1] = sqrt(s[1] / n);
    s[2] = sqrt(s[2] / n);

    return s;
}

__global__ void squareDifferenceKernelUsingByPSNR(float *f0, float *f1, double *tem, int nx, int ny, int nz) {

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x < nx && y < ny && z < nz) {

        tem[arr3(x, y, z)] = f0[arr3(x, y, z)] - f1[arr3(x, y, z)];
        tem[arr3(x, y, z)] *= tem[arr3(x, y, z)];
    }
}

// __global__ void differenceKernelUsingByPSNR(float* f0,float* f1,float* tem,int nx,int ny,int nz){

//     int x=blockIdx.x*blockDim.x+threadIdx.x;
//     int y=blockIdx.y*blockDim.y+threadIdx.y;
//     int z=blockIdx.z*blockDim.z+threadIdx.z;

//     if(x<nx&&y<ny&&z<nz){

//         tem[arr3(x,y,z)]=abs(f0[arr3(x,y,z)]-f1[arr3(x,y,z)]);
//     }
// }

double findMSEn(float *f0, float *f1, int nx, int ny, int nz) {

    int n = nx * ny * nz;

    dim3 blocksz = dim3(BLOCKSZX, BLOCKSZY, BLOCKSZZ);
    dim3 blocknum = dim3(nx / BLOCKSZX + (nx % BLOCKSZX > 0), ny / BLOCKSZY + (ny % BLOCKSZY > 0), nz / BLOCKSZZ + (nz % BLOCKSZZ > 0));

    double *tem;

    cudaMalloc(&tem, n * sizeof(double));

    squareDifferenceKernelUsingByPSNR<<<blocknum, blocksz>>>(f0, f1, tem, nx, ny, nz);

    cudaDeviceSynchronize();

    double ans = sumupUsingReductionForceDouble(tem, n);

    cudaFree(tem);

    return ans;
}

double findPSNR(float *f0, float *f1, int nx, int ny, int nz) {
//double findPSNR(float* f0,float* f1,int nx,int ny,int nz,int output_ae=0){

    int n = nx * ny * nz;

    double Rf0 = findMaximumUsingReduction(f0, n) - findMinimumUsingReduction(f0, n);

    if (Rf0 <= 0) return -inf;

    double MSEn = findMSEn(f0, f1, nx, ny, nz);

    if (MSEn <= 0) return inf;

    double PSNR = 20. * log10(Rf0 * sqrt(n) / sqrt(MSEn));

    // if(output_ae){

    //     float* tem;

    //     cudaMalloc(&tem,n*sizeof(float));

    //     dim3 blocksz=dim3(BLOCKSZX,BLOCKSZY,BLOCKSZZ);
    //     dim3 blocknum=dim3(nx/BLOCKSZX+(nx%BLOCKSZX>0),ny/BLOCKSZY+(ny%BLOCKSZY>0),nz/BLOCKSZZ+(nz%BLOCKSZZ>0));
    //     differenceKernelUsingByPSNR<<<blocknum,blocksz>>>(f0,f1,tem,nx,ny,nz);

    //     float ae=findMaximumUsingReduction(tem,n);

    //     std::cout<<std::scientific<<"Absolute Error="<<ae<<std::endl;
    //     std::cout<<std::fixed;

    //     cudaFree(tem);
    // }

    return PSNR;
}




// device -> device
void derivatives(float *data, int nx, int ny, int nz,
                 float *&dx, float *&dy, float *&dz,
                 float *&dx2, float *&dy2, float *&dz2,
                 float *&dxy, float *&dyz, float *&dzx,
                 float *&gl, float *&lap
) {

    int n = nx * ny * nz;

    cudaMalloc(&dx, n * sizeof(float));
    cudaMalloc(&dy, n * sizeof(float));
    cudaMalloc(&dz, n * sizeof(float));
    cudaMalloc(&dx2, n * sizeof(float));
    cudaMalloc(&dy2, n * sizeof(float));
    cudaMalloc(&dz2, n * sizeof(float));
    cudaMalloc(&dxy, n * sizeof(float));
    cudaMalloc(&dyz, n * sizeof(float));
    cudaMalloc(&dzx, n * sizeof(float));
    cudaMalloc(&gl, n * sizeof(float));
    cudaMalloc(&lap, n * sizeof(float));

    dim3 blocksz = dim3(BLOCKSZX, BLOCKSZY, BLOCKSZZ);
    dim3 blocknum = dim3(nx / BLOCKSZX + (nx % BLOCKSZX > 0), ny / BLOCKSZY + (ny % BLOCKSZY > 0), nz / BLOCKSZZ + (nz % BLOCKSZZ > 0));

    derivativesKernel<<<blocknum, blocksz>>>(data, dx, dy, dz, dx2, dy2, dz2, dxy, dyz, dzx, gl, lap, nx, ny, nz);

    cudaDeviceSynchronize();
}

// host -> device -> host
std::vector<float> derivativesPSNR(float *host_f0, float *host_f1, int nx, int ny, int nz) {

    int n = nx * ny * nz;

    float *f0;
    cudaMalloc(&f0, n * sizeof(float));
    cudaMemcpy(f0, host_f0, n * sizeof(float), cudaMemcpyHostToDevice);

    float *f0_dx, *f0_dy, *f0_dz;
    float *f0_dx2, *f0_dy2, *f0_dz2;
    float *f0_dxy, *f0_dyz, *f0_dzx;
    float *f0_gl, *f0_lap;

    derivatives(f0, nx, ny, nz, f0_dx, f0_dy, f0_dz, f0_dx2, f0_dy2, f0_dz2, f0_dxy, f0_dyz, f0_dzx, f0_gl, f0_lap);
    std::vector<float> f0_sobolev = sobolev(f0, f0_dx, f0_dy, f0_dz, f0_dx2, f0_dy2, f0_dz2, f0_dxy, f0_dyz, f0_dzx, nx, ny, nz);

    float *f1;
    cudaMalloc(&f1, n * sizeof(float));
    cudaMemcpy(f1, host_f1, n * sizeof(float), cudaMemcpyHostToDevice);

    float *f1_dx, *f1_dy, *f1_dz;
    float *f1_dx2, *f1_dy2, *f1_dz2;
    float *f1_dxy, *f1_dyz, *f1_dzx;
    float *f1_gl, *f1_lap;

    derivatives(f1, nx, ny, nz, f1_dx, f1_dy, f1_dz, f1_dx2, f1_dy2, f1_dz2, f1_dxy, f1_dyz, f1_dzx, f1_gl, f1_lap);
    std::vector<float> f1_sobolev = sobolev(f1, f1_dx, f1_dy, f1_dz, f1_dx2, f1_dy2, f1_dz2, f1_dxy, f1_dyz, f1_dzx, nx, ny, nz);

    std::vector<float> vec(12);

    vec[0] = findPSNR(f0, f1, nx, ny, nz);

    vec[1] = findPSNR(f0_dx, f1_dx, nx, ny, nz);
    vec[2] = findPSNR(f0_dy, f1_dy, nx, ny, nz);
    vec[3] = findPSNR(f0_dz, f1_dz, nx, ny, nz);

    vec[4] = findPSNR(f0_dx2, f1_dx2, nx, ny, nz);
    vec[5] = findPSNR(f0_dy2, f1_dy2, nx, ny, nz);
    vec[6] = findPSNR(f0_dz2, f1_dz2, nx, ny, nz);

    vec[7] = findPSNR(f0_gl, f1_gl, nx, ny, nz);
    vec[8] = findPSNR(f0_lap, f1_lap, nx, ny, nz);

    vec[9] = std::abs(f0_sobolev[0] - f1_sobolev[0]) / f0_sobolev[0];
    vec[10] = std::abs(f0_sobolev[1] - f1_sobolev[1]) / f0_sobolev[1];
    vec[11] = std::abs(f0_sobolev[2] - f1_sobolev[2]) / f0_sobolev[2];

    cudaFree(f0);
    cudaFree(f0_dx), cudaFree(f0_dy), cudaFree(f0_dz);
    cudaFree(f0_dx2), cudaFree(f0_dy2), cudaFree(f0_dz2);
    cudaFree(f0_dxy), cudaFree(f0_dyz), cudaFree(f0_dzx);
    cudaFree(f0_gl), cudaFree(f0_lap);

    cudaFree(f1);
    cudaFree(f1_dx), cudaFree(f1_dy), cudaFree(f1_dz);
    cudaFree(f1_dx2), cudaFree(f1_dy2), cudaFree(f1_dz2);
    cudaFree(f1_dxy), cudaFree(f1_dyz), cudaFree(f1_dzx);
    cudaFree(f1_gl), cudaFree(f1_lap);

    return vec;
}

#undef BLOCKSZX
#undef BLOCKSZY
#undef BLOCKSZZ
#undef blksz