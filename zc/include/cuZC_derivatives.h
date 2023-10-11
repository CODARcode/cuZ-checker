#ifndef _CUZC_DERIVATIVES_H
#define _CUZC_DERIVATIVES_H
#include <vector>
/**
 * Compute the derivatives for 3D data
 * Derivatives are the same size as the inpout
 * @param data 3D data
 * @param nx data shape: data[dx][dy][dz]
 * @param ny Data shape: data[dx][dy][dz]
 * @param nz Data shape: data[dx][dy][dz]
 * @param dx output, partial derivative dx
 * @param dy output, partial derivative dy
 * @param dz output, partial derivative dz
 * @param dx2 output, second order derivative dx^2
 * @param dy2 output, second order derivative dy^2
 * @param dz2 output, second order derivative dz^2
 * @param dxy output, second order derivative dxdy
 * @param dyz output, second order derivative dxdz
 * @param dzx output, second order derivative dxdz
 * @param gl output, gradient length
 * @param lap output, laplacian
 */
void derivatives(float *data, int nx, int ny, int nz,
                 float *&dx, float *&dy, float *&dz,
                 float *&dx2, float *&dy2, float *&dz2,
                 float *&dxy, float *&dyz, float *&dzx,
                 float *&gl, float *&lap
);

/**
 * Compute the derivatives for two 3D data (data1, data2), and output the PSNR of the derivatives of the two
 * @param data1 3D data 1
 * @param data2 3D data 2
 * @param nx data shape: data1[dx][dy][dz] and data2[dx][dy][dz]
 * @param ny data shape: data1[dx][dy][dz] and data2[dx][dy][dz]
 * @param nz data shape: data1[dx][dy][dz] and data2[dx][dy][dz]
 * @return a vector contains: (PSNR of data, PSNR of dx,dy,dz, PSNR of dx2,dy2,dz2, PSNR of gradient length, laplacian, relative error of s0,s1,s2)
 */
std::vector<float> derivativesPSNR(float *data1, float *data2, int nx, int ny, int nz);

#endif