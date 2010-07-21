#include <cuda.h>
#include <math.h>
#include "stdio.h"
#include "gpu_kernel.h"
#include "assert.h"
#include "c_common.h"

char buf[100000];

__global__ void GPUKernel(float *Lm, float *Um, float *LUm,
                          int *Lstart, int *Lend,
                          int *Ustart, int *Uend,
                          int block, int total) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int lstart = Lstart[i], ustart = Ustart[i];
  float val = 0.0;
  int offsetx = (i % (block * block)) % block;
  int offsety = (i % (block * block)) / block;

  for (int k = 0; k < block; k++) {
    int i1 = (k + block * offsety) + lstart;
    int i2 = (k * block + offsetx) + ustart;

    float elm1 = Lm[i1];
    float elm2 = Um[i2];
    val += elm1 * elm2;
  }

  LUm[i] += -1 * val;
}

void checkCUDAError(const char *msg) {
  cudaError_t err = cudaGetLastError();
  if (cudaSuccess != err) {
    fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
    exit(EXIT_FAILURE);
  }            
}

void checkCUDAKernelError(const char *msg) {
  cudaError_t err = cudaThreadSynchronize();
  if( cudaSuccess != err) {
    sprintf(buf, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
    fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
    exit(EXIT_FAILURE);
  }                        
}

static int n = MAX_OFFLOAD_SIZE;
static float *d_Lm;
static float *d_Um;
static float *d_LUm;
static int *d_Lstart;
static int *d_Lend;
static int *d_Ustart;
static int *d_Uend;
static int allocated = 0;
  
void GPUallocate() {
  if (!allocated) {
    printf("allocation happening\n");
    cudaMalloc((void **) &d_Lm, n * sizeof(float));
    cudaMalloc((void **) &d_Um, n * sizeof(float));
    cudaMalloc((void **) &d_LUm, n * sizeof(float));
    cudaMalloc((void **) &d_Lstart, n * sizeof(int));
    cudaMalloc((void **) &d_Lend, n * sizeof(int));
    cudaMalloc((void **) &d_Ustart, n * sizeof(int));
    cudaMalloc((void **) &d_Uend, n * sizeof(int));
    checkCUDAError("mem allocation");
    allocated = 1;
    printf("allocation finished\n");
  }
}

void GPUKernelDGEMM(float Lm[], float Um[], float LUm[], 
                    int Lstart[], int Lend[], int Ustart[], int Uend[],
                    int block, int total) {
  /*for (int i = 0; i < sn; i++) {
    assert(sStart[i] >= 0);
    assert(sStart[i] <= sEnd[i]);
    assert(sEnd[i] < fn);
    assert(sEnd[i] - sStart[i] > 0 && sEnd[i] - sStart[i] < 5000);
  }

  for (int i = 0; i < fn; i++) {
    assert(fStart[i] >= 0);
    assert(fStart[i] <= fEnd[i]);
    assert(fEnd[i] < sn);
    assert(fEnd[i] - fStart[i] > 0 && fEnd[i] - fStart[i] < 5000);
    }*/
	
  /*int findex = fn - 1;
    int sindex = sn - 1;*/

  //printf("fStart[0] = %d, fEnd[0] = %d, fStart[%d] = %d, fEnd[%d] = %d\n", fStart[0], fEnd[0], findex, fStart[findex], findex, fEnd[findex]);
  //printf("sStart[0] = %d, sEnd[0] = %d, sStart[%d] = %d, sEnd[%d] = %d\n", sStart[0], sEnd[0], sindex, sStart[sindex], sindex, sEnd[sindex]);

  //printf("total = %d, block = %d\n", total, block);

  size_t dsize = total * sizeof(float);
  size_t isize = total * sizeof(int);

  //printf("dsize = %d, isize = %d\n", (int)dsize, (int)isize);
	
  /*assert(sizes < (n * sizeof(float)) );
    assert(sizef < (n * sizeof(float)) );
    assert(sizesi < (n * sizeof(int)) );
    assert(sizefi < (n * sizeof(int)) );*/

  GPUallocate();

  cudaMemcpy(d_Lm, Lm, dsize, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Um, Um, dsize, cudaMemcpyHostToDevice);
  cudaMemcpy(d_LUm, LUm, dsize, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Lstart, Lstart, isize, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Lend, Lend, isize, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Ustart, Ustart, isize, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Uend, Uend, isize, cudaMemcpyHostToDevice);

  checkCUDAError("mem copy to device");

  /*for (int i = 0; i < total; i++) {
    printf("before LUm[%d] = %f\n", i, LUm[i]);
  }

  for (int i = 0; i < total; i++) {
    printf("Lstart[%d] = %d\n", i, Lstart[i]);
  }

  for (int i = 0; i < total; i++) {
    printf("Ustart[%d] = %d\n", i, Ustart[i]);
  }

  for (int i = 0; i < total; i++) {
    printf("before Lm[%d] = %f\n", i, Lm[i]);
  }

  for (int i = 0; i < total; i++) {
    printf("before Um[%d] = %f\n", i, Um[i]);
  }

  printf("TOTAL = %d\n", total);*/

  int blockSize = 16;
  int nBlocks = total/blockSize + (total%blockSize == 0?0:1);

  GPUKernel<<<nBlocks, blockSize>>>(d_Lm, d_Um, d_LUm,
                                    d_Lstart, d_Lend,
                                    d_Ustart, d_Uend,
                                    block, total);
  checkCUDAKernelError("kernel execute");

  cudaMemcpy(LUm, d_LUm, dsize, cudaMemcpyDeviceToHost);

  /*for (int i = 0; i < total; i++) {
    printf("after LUm[%d] = %f\n", i, LUm[i]);
    }*/

  checkCUDAError("mem copy from device");
}
