#include<iostream>
#include<iomanip>
#include<fstream>
#include<assert.h>
#include"cuZC_derivatives.h"

using namespace std;

template<typename Type>
void readFile(const char *file, const size_t num, Type *data) {
    std::ifstream fin(file, std::ios::binary);
    if (!fin) {
        std::cout << " Error, Couldn't find the file: " << file << "\n";
        exit(0);
    }
    fin.seekg(0, std::ios::end);
    const size_t num_elements = fin.tellg() / sizeof(Type);
    assert(num_elements == num && "File size is not equals to the input setting");
    fin.seekg(0, std::ios::beg);
    fin.read(reinterpret_cast<char *>(data), num_elements * sizeof(Type));
    fin.close();
}

template<typename Type>
void writeFile(const char *file, const size_t num_elements, Type *data) {
    std::ofstream fout(file, std::ios::binary);
    fout.write(reinterpret_cast<const char *>(&data[0]), num_elements * sizeof(Type));
    fout.close();
}

char f0_file_name[64], f1_file_name[64], buf[64];
int nx, ny, nz, n;

float *f0, *f1;

signed main() {

    printf("please input the original file name:\n");

    scanf("%s", f0_file_name);

    printf("please input the unzipped file name:\n");

    scanf("%s", f1_file_name);

    printf("please input the size like [nx] [ny] [nz]:\n");

    scanf("%d %d %d", &nx, &ny, &nz);

    n = nx * ny * nz;

    f0 = (float *) malloc(n * sizeof(float));
    f1 = (float *) malloc(n * sizeof(float));

    readFile<float>(f0_file_name, n, f0);
    readFile<float>(f1_file_name, n, f1);

    vector<float> vec = derivativesPSNR(f0, f1, nx, ny, nz);

    cout << "PSNR of data = " << vec[0] << endl;

    cout << "PSNR of dx = " << vec[1] << endl;
    cout << "PSNR of dy = " << vec[2] << endl;
    cout << "PSNR of dz = " << vec[3] << endl;

    cout << "PSNR of dx2 = " << vec[4] << endl;
    cout << "PSNR of dy2 = " << vec[5] << endl;
    cout << "PSNR of dz2 = " << vec[6] << endl;

    cout << "PSNR of gradient length = " << vec[7] << endl;
    cout << "PSNR of laplacian = " << vec[8] << endl;

    cout << scientific << setprecision(2);
    cout << "Relative Error of s0 = " << vec[9] << endl;
    cout << "Relative Error of s1 = " << vec[10] << endl;
    cout << "Relative Error of s2 = " << vec[11] << endl;

    free(f0), free(f1);

    return 0;
}