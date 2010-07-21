#ifndef GPU_OFFLOADER
#define GPU_OFFLOADER

#include <vector>
#include "gpuwork.decl.h"

using namespace std;

class GPUWork : public CBase_GPUWork {
public:
  GPUWork() {}
  GPUWork(CkMigrateMessage* msg) {}

  void gpu_offload(list<JMessage*>& msgs);
  void FakeGPUDGEMM(int tsize, int agglom, int block, float *Ap, float *Bp,
                    float *Cp, int *Astart, int *Bstart);
};

#endif
