#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <fstream>
#include <iostream>
#include <cassert>
#include <cinttypes>
#include "cuZC_entry.h"

using namespace std;

template<typename Type>
void readFile(const char *file, const size_t num, Type *data) {
    ifstream fin(file, ios::binary);
    if (!fin) {
        cout << "Error, Couldn't find the file: " << file << "\n";
        exit(1);
    }
    fin.seekg(0, ios::end);
    const size_t num_elements = fin.tellg() / sizeof(Type);
    assert(num_elements == num && "File size does not match the input dimensions");
    fin.seekg(0, ios::beg);
    fin.read(reinterpret_cast<char *>(data), num_elements * sizeof(Type));
    fin.close();
}

double cpu_ssim_3d(float *data1, float *data2,
                   int r3, int r2, int r1,
                   int ssimSize, int ssimShift)
{
    int nw = 0;
    double ssimSum = 0.0;

    for (int oz = 0; oz + ssimSize <= r1; oz += ssimShift) {
        for (int ox = 0; ox + ssimSize <= r3; ox += ssimShift) {
            for (int oy = 0; oy + ssimSize <= r2; oy += ssimShift) {
                float xMin = FLT_MAX, xMax = -FLT_MAX;
                float yMin = FLT_MAX, yMax = -FLT_MAX;
                double xSum = 0, x2Sum = 0;
                double ySum = 0, y2Sum = 0;
                double xySum = 0;
                int np = 0;

                for (int wx = 0; wx < ssimSize; wx++) {
                    for (int wy = 0; wy < ssimSize; wy++) {
                        for (int wz = 0; wz < ssimSize; wz++) {
                            int idx = (ox + wx) * r1 * r2
                                    + (oy + wy) * r1
                                    + (oz + wz);
                            float x = data1[idx];
                            float y = data2[idx];
                            np++;
                            if (x < xMin) xMin = x;
                            if (x > xMax) xMax = x;
                            if (y < yMin) yMin = y;
                            if (y > yMax) yMax = y;
                            xSum  += x;
                            x2Sum += (double)x * x;
                            ySum  += y;
                            y2Sum += (double)y * y;
                            xySum += (double)x * y;
                        }
                    }
                }

                double xMean  = xSum / np;
                double yMean  = ySum / np;
                double xSigma = sqrt(fabs(x2Sum / np - xMean * xMean));
                double ySigma = sqrt(fabs(y2Sum / np - yMean * yMean));
                double xyCov  = xySum / np - xMean * yMean;

                double c1, c2;
                if (xMax - xMin == 0) {
                    c1 = 0.01 * 0.01;
                    c2 = 0.03 * 0.03;
                } else {
                    double range = xMax - xMin;
                    c1 = 0.01 * 0.01 * range * range;
                    c2 = 0.03 * 0.03 * range * range;
                }
                double c3 = c2 / 2;

                double lum  = (2 * xMean * yMean + c1) / (xMean * xMean + yMean * yMean + c1);
                double cont = (2 * xSigma * ySigma + c2) / (xSigma * xSigma + ySigma * ySigma + c2);
                double stru = (xyCov + c3) / (xSigma * ySigma + c3);
                ssimSum += lum * cont * stru;
                nw++;
            }
        }
    }

    printf("CPU: %d windows, SSIM sum = %e\n", nw, ssimSum);
    return ssimSum / nw;
}

int main(int argc, char *argv[])
{
    if (argc < 6) {
        printf("Usage: test_cussim ori_file dec_file r3 r2 r1\n");
        printf("  Dimensions are slowest-to-fastest (r3 x r2 x r1).\n");
        printf("  Example: test_cussim orig.dat dec.dat 128 8 8\n");
        printf("  Use the same file twice for a sanity check (SSIM should be 1.0).\n");
        return 1;
    }

    const char *oriFile = argv[1];
    const char *decFile = argv[2];
    size_t r3 = strtoimax(argv[3], NULL, 10);
    size_t r2 = strtoimax(argv[4], NULL, 10);
    size_t r1 = strtoimax(argv[5], NULL, 10);
    size_t n  = r3 * r2 * r1;

    int ssimSize  = 7;
    int ssimShift = 1;
    if (argc >= 7) ssimSize  = atoi(argv[6]);
    if (argc >= 8) ssimShift = atoi(argv[7]);

    printf("ori: %s\n", oriFile);
    printf("dec: %s\n", decFile);
    printf("dimensions: r3=%zu  r2=%zu  r1=%zu  (total %zu elements)\n", r3, r2, r1, n);
    printf("SSIM window: size=%d  shift=%d\n", ssimSize, ssimShift);

    if ((int)r1 < ssimSize || (int)r2 < ssimSize || (int)r3 < ssimSize) {
        printf("Error: all dimensions must be >= ssimSize (%d)\n", ssimSize);
        return 1;
    }

    float *data1 = (float *)malloc(n * sizeof(float));
    float *data2 = (float *)malloc(n * sizeof(float));
    readFile<float>(oriFile, n, data1);
    readFile<float>(decFile, n, data2);

    // --- GPU SSIM ---
    printf("\n=== GPU SSIM ===\n");
    cu_SSIM(data1, data2, r3, r2, r1, ssimSize, ssimShift);

    int blksize = (r1 - ssimSize) / ssimShift + 1;
    int nw_r3   = (r3 - ssimSize) / ssimShift + 1;
    int nw_r2   = (r2 - ssimSize) / ssimShift + 1;
    int totalWindows = blksize * nw_r3 * nw_r2;
    printf("(total windows = %d)\n", totalWindows);

    // --- CPU reference SSIM ---
    printf("\n=== CPU SSIM (reference) ===\n");
    double cpuSSIM = cpu_ssim_3d(data1, data2, r3, r2, r1, ssimSize, ssimShift);
    printf("CPU average SSIM = %.10f\n", cpuSSIM);

    free(data1);
    free(data2);
    return 0;
}
