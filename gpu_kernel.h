#ifndef GPU_KERNEL_H
#define GPU_KERNEL_H

const float EPS = 0.1;

#define ASYNC 0

void GPUallocate();
void GPUSingleallocate();

void GPUKernelDGEMM(float Lm[], float Um[], float LUm[], 
                    int Lstart[], int Lend[], int Ustart[], int Uend[],
                    int block, int total);
void GPUSingleOffload(float Lm[], float Um[], float LUm[], int block);

#endif
