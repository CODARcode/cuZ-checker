#ifndef CUZC_ENTRY_H
#define CUZC_ENTRY_H

#include <stdio.h>
#include "matrix.hpp"
#include "cuZC_ssim.h"
#include "cuZC_typeOne.h"
#include "cuZC_typeTwo.h"
#include "cuZC_typeThree.h"

double cu_SSIM_3d_windowed(int windowSize0, int windowSize1, int windowSize2, int windowShift0, int windowShift1, int windowShift2);
int cu_SSIM(float *data1, float *data2, size_t r3, size_t r2, size_t r1, int ssimSize, int ssimShift);

double* cu_typeOne(float *ddata1, float *ddata2, double *ddiff, double *absErrPDF, double *hresults, size_t r3, size_t r2, size_t r1, size_t ne);
float *cu_typeTwo(float *ddata, float *der, size_t r3, size_t r2, size_t r1, double avg, size_t order);
double cu_typeThree(float *data1, float *data2, int r3, int r2, int r1, int ssimSize, int ssimShift); 

#endif /* ----- #ifndef CUZC_ENTRY_H  ----- */
