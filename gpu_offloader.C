#include "pup_stl.h"
#include "c_common.h"
#include "Scheduler.h"
#include "gpu_kernel.h"
#include "assert.h"

#include "gpuwork.decl.h"

using namespace std;

extern CProxy_Scheduler scheduler;

class GPUWork : public CBase_GPUWork {
public:
  GPUWork() {}
  GPUWork(CkMigrateMessage* msg) {}

  // L same row, U same col
  void gpu_offload(list<JMessage>& msgs) {
    //ckout << "---- offloading to GPU" << endl;

    /*char note[200];
      sprintf(note, "gpu offloading happening");
      traceUserSuppliedNote(note);*/
    
    int asize = msgs.front().fsize * msgs.front().fsize;

    //CkPrintf("asize = %d\n", asize);
      
    float *Lvec = new float[msgs.size() * asize];
    float *Uvec= new float[msgs.size() * asize];
    float *Avec= new float[msgs.size() * asize];

    int *Lstart = new int[msgs.size() * asize];
    int *Lend = new int[msgs.size() * asize];
    int *Ustart = new int[msgs.size() * asize];
    int *Uend = new int[msgs.size() * asize];

    int Li = 0, Ui = 0, Ai = 0, 
      Lsi = 0, Lei = 0, Usi = 0, Uei = 0;

    int lastX = -1, lastY = -1;
    int startL = 0, startU = 0;

    int block = 0;

    for (list<JMessage>::iterator iter = msgs.begin(); 
	 iter != msgs.end(); ++iter) {
      bool copyL = true, copyU = true;
      int tsize = iter->fsize * iter->ssize;

      block = iter->fsize;

      //if (lastX == iter->x)  
      //copyL = false;
      if (lastX != -1)
        startL += tsize;

      //if (lastY == iter->y)
      //copyU = false;
      if (lastY != -1)
         startU += tsize;

      for (int i = 0; i < tsize; i++) {
        if (copyL)
          Lvec[Li++] = iter->first[i];

        if (copyU)
          Uvec[Ui++] = iter->second[i];

        Avec[Ai++] = iter->LU[i];

        Lstart[Lsi++] = startL;
        Lend[Lei++] = startL + tsize;

        Ustart[Usi++] = startU;
        Uend[Uei++] = startU + tsize;
      }

      lastX = iter->x;
      lastY = iter->y;
    }

    int size = Ui;

    //CkPrintf("AA size = %d\n", size);

    //float *AmP = (float*)malloc(sizeof(float) * size);
    //memcpy(AmP, &Avec[0], sizeof(float) * size);

    GPUKernelDGEMM(Lvec, Uvec, Avec, Lstart, Ustart, block, size);

    /*FakeGPUDGEMM(size, msgs.size(), block, Uvec, Lvec, Avec,
      Ustart, Uend, Lstart, Lend);*/

    /*for (int i = 0; i < size; i++) {
      if (Am[i] != AmP[i])
        ckout << "FOUND: " << i << " Am = " << Am[i] << ", AmP = " << AmP[i] <<
          ", Lstart = " << Ls[i] << endl;
          }*/

    int firstLoc = 0;

    for (list<JMessage>::iterator iter = msgs.begin();
	 iter != msgs.end(); ++iter) {
      int tsize = iter->fsize * iter->ssize;

      for (int i = 0; i < tsize; i++) {
        iter->LU[i] = Avec[i + firstLoc];
      }

      firstLoc += tsize;
    }

    delete[] Lvec, Uvec, Avec, Lstart, Lend, Ustart, Uend;

    scheduler[CkMyPe()].finishedGPU(msgs);
  }

  void FakeGPUDGEMM(int tsize, int agglom, int block, float *Ap, float *Bp,
                    float *Cp, int *Astart, int *Aend, int *Bstart,
                    int *Bend) {
    int bsize = tsize / agglom;

    for (int m = 0; m < agglom; m++) {
      /*ckout << "bstart = " << Bstart[bsize * m] << endl;
        ckout << "astart = " << Astart[bsize * m] << endl;*/

      float *A = Astart[bsize * m] + Ap;
      float *B = Bstart[bsize * m] + Bp;
      float *C = m * bsize + Cp;

      /*for (int i = 0; i < block; i++) {
        ckout << "L = " << B[i] << endl;
        ckout << "U = " << A[i] << endl;
        ckout << "LU = " << C[i] << endl;
        }*/

      int N = block;
      unsigned i, j, k;

      for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) {
          const float *Ai_ = A + i;
          const float *B_j = B + j*N;

          float cij = *(C + j*N + i);

          for (k = 0; k < N; ++k) {
            cij += -1 * *(Ai_ + k*N) * *(B_j + k);
          }

          *(C + j*N + i) = cij;
        }
      }
    }
  }
};

#include "gpuwork.def.h"
