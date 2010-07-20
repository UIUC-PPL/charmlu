#ifndef GPU_KERNEL_H
#define GPU_KERNEL_H

const float EPS = 0.1;

#define ASYNC 0

void GPUallocate();

void GPUKernelDGEMM(float Lm[], float Um[], float LUm[], 
                    int Lstart[], int Lend[], int Ustart[], int Uend[],
                    int block, int total);

#endif
