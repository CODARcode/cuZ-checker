#ifndef CUZC_ENTRY_H
#define CUZC_ENTRY_H

#include <stdio.h>
#include "matrix.hpp"
#include "cuZC_ssim.h"
 
double cu_SSIM_3d_windowed(int windowSize0, int windowSize1, int windowSize2, int windowShift0, int windowShift1, int windowShift2);
int cu_SSIM(float *data1, float *data2, size_t r3, size_t r2, size_t r1, int ssimSize, int ssimShift);

#endif /* ----- #ifndef CUZC_ENTRY_H  ----- */
