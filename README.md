# cuZ-checker
This is the cuda GPU version of Z-checker (https://github.com/CODARcode/Z-checker), a library to characterize the data and check the compression results of lossy compressors

## Prerequisites
cmake >= 3.19 <br />
CUDA toolkit >= 10.0 with samples/common/inc

## How to install
Under the cuZ-checker directory:
```bash
mkdir build
cd build
cmake ../
make -j8
```

## How to run
Assuming your original 3D (100\*500\*500) data is foo.dat, the compressed data is foo.dat.sz, and the decompressed data is foo.dat.sz.out. They all are under cuZ-checker/example/Datasets. <br />
To perform the measurements, enter cuZ-checker/example and first run data property analysis:
```bash
./analyzeDataProperty foo -f zc.config Datasets/foo.dat 500 500 100
```
Then conduct the comperisons between the original and decompressed data:
```bash
./compareDataSets -f zc.config "SZ(1E-3)" foo Datasets/foo.dat Datasets/foo.dat.sz.out 500 500 100
```
Note that "SZ(1E-3)" is the error bound that should be consistent with the compressor. <br />
zc.config is changable that controls the measurement settings.
