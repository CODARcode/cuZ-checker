#include "cuZC_entry.h"
#include "timingGPU.h"

TimingGPU timer_GPU;
 
double cu_SSIM_3d_windowed(int windowSize0, int windowSize1, int windowSize2, int windowShift0, int windowShift1, int windowShift2)
{
    char a[N] = "Hello \0\0\0\0\0\0";
    int b[N] = {15, 10, 6, 0, -11, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
     
    char *ad;
    int *bd;
    const int csize = N*sizeof(char);
    const int isize = N*sizeof(int);
     
    printf("%s", a);
     
    cudaMalloc( (void**)&ad, csize  ); 
    cudaMalloc( (void**)&bd, isize  ); 
    cudaMemcpy( ad, a, csize, cudaMemcpyHostToDevice  ); 
    cudaMemcpy( bd, b, isize, cudaMemcpyHostToDevice  ); 

    dim3 dimBlock( blocksize, 1  );
    dim3 dimGrid( 1, 1  );
    //hello<<<dimGrid, dimBlock>>>(ad, bd);
    cudaMemcpy( a, ad, csize, cudaMemcpyDeviceToHost  ); 
    cudaFree( ad  );
    cudaFree( bd  );

    printf("%s\n", a);
}

int cu_SSIM(float *data1, float *data2, size_t r3, size_t r2, size_t r1, int ssimSize, int ssimShift)
{
    float data[246];
    for (int i=0; i<246; i++){
        data[i] = 1;
    }
    int blksize = (r1 - ssimSize) / ssimShift + 1;
    int xsize = ((r2 - ssimSize) / ssimShift + 1)*((r3 - ssimSize) / ssimShift + 1);

    double results[blksize] = { 0 };
    //printf("test=%f, %f\n", data[32], results[32]);

    float *ddata1, *ddata2;
    double *dresults;
    //for (int i=r1*r2*6+r2*6;i<r1*r2*6+r2*6+7;i++){
    ////for (int i=0;i<r1*r2*r3;i++){
    //    printf("data%i=%e, %e\n",i, data1[i], data2[i]);
    //    printf("data%i=%e, %e\n",i, data1[i], data2[i]);

    //}

    const int csize = r3 * r2 * r1 * sizeof(float);
    const int isize = blksize * sizeof(double);

    cudaMalloc((void**)&ddata1,   csize); 
    cudaMalloc((void**)&ddata2,   csize); 
    cudaMalloc((void**)&dresults, isize); 
    cudaMemcpy(ddata1,   data1,   csize, cudaMemcpyHostToDevice); 
    cudaMemcpy(ddata2,   data2,   csize, cudaMemcpyHostToDevice); 
    cudaMemcpy(dresults, results, isize, cudaMemcpyHostToDevice); 

    timer_GPU.StartCounter();
    dim3 dimBlock(64, 1);
    dim3 dimGrid(blksize, 1);
    ssim<<<dimGrid, dimBlock>>>(ddata1, ddata2, dresults, r3, r2, r1, ssimSize, ssimShift);
    cudaMemcpy(results, dresults, isize, cudaMemcpyDeviceToHost); 
    double x=0;
    printf("GPU timing: %f ms\n", timer_GPU.GetCounter());
    for (int i=0; i<blksize; i++){
        x += results[i];
        printf("results%i=%e\n",i,x);

    }

    cudaFree(ddata1);
    cudaFree(ddata2);
    cudaFree(dresults);

    return 0;
}

double *cu_typeOne(float *ddata1, float *ddata2, double *ddiff, double *absErrPDF, double *results, size_t r3, size_t r2, size_t r1, size_t ne){

    //float *ddata1, *ddata2;
    double *dabsErrPDF, *dresults;
    //for (int i=r1*r2*6+r2*6;i<r1*r2*6+r2*6+7;i++){
    ////for (int i=0;i<r1*r2*r3;i++){
    //    printf("data%i=%e, %e\n",i, data1[i], data2[i]);
    //    printf("data%i=%e, %e\n",i, data1[i], data2[i]);

    //}

    const int dsize = ne * sizeof(double);
    const int rsize = r3 * 10 * sizeof(double);

    cudaMalloc((void**)&dabsErrPDF, dsize); 
    cudaMalloc((void**)&dresults, rsize); 
    cudaMemcpy(dresults, results, rsize, cudaMemcpyHostToDevice); 

    timer_GPU.StartCounter();
    void *kernelArgs[] = {
        (void *)&ddata1, (void *)&ddata2, (void *)&ddiff, (void *)&dresults, 
        (void *)&r3, (void *)&r2, (void *)&r1, (void *)&ne,
    };
    dim3 dimBlock(32, 10);
    dim3 dimGrid(r3, 1);
    const int sMemsize = 10 * dimBlock.y * sizeof(double);
    printf("test:%i,%i\n",sMemsize, dimBlock.y);
    //type_one<<<dimGrid, dimBlock>>>(ddata1, ddata2, ddiff, dresults, r3, r2, r1, ne);
    checkCudaErrors(cudaLaunchCooperativeKernel((void*)type_one,
                                                dimGrid, dimBlock, kernelArgs, sMemsize));

    //dim3 dimBlock2(32, 10);
    //gridReduction<<<1, dimBlock2>>>(dresults, r3);

    cudaMemcpy(results, dresults, rsize, cudaMemcpyDeviceToHost); 
    double x=0;
    printf("GPU timing: %f ms\n", timer_GPU.GetCounter());
    //for (int i=0; i<r3; i++){
    //    x += results[i];
    //    printf("results%i=%e\n",i,x);

    //}

    cudaFree(dabsErrPDF);
    cudaFree(dresults);

    return results;
}

float *cu_typeTwo(float *ddata, float *der, size_t r3, size_t r2, size_t r1, double avg, size_t order){

    float *dder, *autocor, *dautocor;
    const int dsize = (r3-order*2) * (r2-order*2) * (r1-order*2) * sizeof(float);
    int gridsize = (r3-order*2)/(16-order*2)+((r3-order*2)%(16-order*2)?1:0);
    int corsize = gridsize * order * 2 * sizeof(float);

    cudaMalloc((void**)&dder, dsize); 
    cudaMemcpy(dder, der, dsize, cudaMemcpyHostToDevice); 
    autocor = (float*)malloc(corsize);
    memset(autocor, 0, corsize);
    cudaMalloc((void**)&dautocor, corsize); 
    cudaMemcpy(dautocor, autocor, corsize, cudaMemcpyHostToDevice); 

    timer_GPU.StartCounter();
    dim3 dimBlock(16, 16);
    dim3 dimGrid(gridsize, 1);
    const int sMemsize = (16 * dimBlock.x * dimBlock.y + dimBlock.y * order * 2) * sizeof(double);
    //type_two<<<dimGrid, dimBlock>>>(ddata, dder, r3, r2, r1, order);
    void *kernelArgs[] = {
        (void *)&ddata, (void *)&dder, (void *)&dautocor,
        (void *)&r3, (void *)&r2, (void *)&r1, (void *)&avg, (void *)&order,
    };
    checkCudaErrors(cudaLaunchCooperativeKernel((void*)type_two,
                                                dimGrid, dimBlock, kernelArgs, sMemsize));

    cudaMemcpy(der, dder, dsize, cudaMemcpyDeviceToHost); 
    cudaMemcpy(autocor, dautocor, corsize, cudaMemcpyDeviceToHost); 

    printf("GPU timing: %f ms\n", timer_GPU.GetCounter());
    //for (int i=0;i<(r3-4)*(r2-4)*(r1-4);i++){
    //    if (der[i]!=0.0) printf("ddata%i=%e\n",i,der[i]);
    //}

    cudaFree(dder);
    cudaFree(dautocor);

    return autocor;
}

double cu_typeThree(float *data1, float *data2, int r3, int r2, int r1, int ssimSize, int ssimShift) 
{

    //for (int i=r1*r2*6+r2*6;i<r1*r2*6+r2*6+7;i++){
    ////for (int i=0;i<r1*r2*r3;i++){
    //    printf("data%i=%e, %e\n",i, data1[i], data2[i]);
    //    printf("data%i=%e, %e\n",i, data1[i], data2[i]);

    //}
    float *ddata1, *ddata2;
    double *dresults;
	double *results = (double*)malloc(r2*r1*((r3-ssimSize+1)/ssimShift)*sizeof(double));
	//double *results = (double*)malloc(r2*r1*sizeof(double));
    memset(results, 0, r2*r1*((r3-ssimSize+1)/ssimShift)*sizeof(double));
    //for (int i=r1*r2*6+r2*6;i<r1*r2*6+r2*6+7;i++){
    ////for (int i=0;i<r1*r2*r3;i++){
    //    printf("data%i=%e, %e\n",i, data1[i], data2[i]);
    //    printf("data%i=%e, %e\n",i, data1[i], data2[i]);

    //}
    int yNum = (8-ssimSize)/ssimShift+1;
    int yTotal = (r2-ssimSize)/ssimShift+1;
    int blksize = yTotal/yNum + (yTotal%yNum?1:0);

    const int csize = r3 * r2 * r1 * sizeof(float);

    cudaMalloc((void**)&ddata1,   csize); 
    cudaMalloc((void**)&ddata2,   csize); 
    cudaMemcpy(ddata1,   data1,   csize, cudaMemcpyHostToDevice); 
    cudaMemcpy(ddata2,   data2,   csize, cudaMemcpyHostToDevice); 

    const int rsize = r1 * r2*((r3-ssimSize+1)/ssimShift) * sizeof(double);

    cudaMalloc((void**)&dresults, rsize); 
    cudaMemcpy(dresults, results, rsize, cudaMemcpyHostToDevice); 

    timer_GPU.StartCounter();
    dim3 dimBlock(32, 8);
    dim3 dimGrid(blksize, 1);
    type_three<<<dimGrid, dimBlock>>>(ddata1, ddata2, dresults, r3, r2, r1, ssimSize, ssimShift, yNum);

    cudaMemcpy(results, dresults, rsize, cudaMemcpyDeviceToHost); 
    double x=0, y=0;
    printf("GPU timing: %f ms\n", timer_GPU.GetCounter());
    for (int i=0; i<r1*r2*((r3-ssimSize+1)/ssimShift); i++){
        if (i%r1==0) {
printf("results%i=%e\n",i/r1-1,x);
x=0;

        }
        //if (i%r1==0) printf("delimiter%i\n",(i/r1));
        x += results[i];
        y += results[i];
        //printf("results%i=%e\n",i,results[i]);

    }
    printf("results=%e\n",y);

    cudaFree(dresults);

    return 0;
}
