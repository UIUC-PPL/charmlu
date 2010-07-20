#include "pup_stl.h"
#include "c_common.h"
#include "Scheduler.h"
#include "assert.h"

#include "gpuwork.decl.h"

using namespace std;

static double A = 2.0; // Force Calculation parameter 1
static double B = 1.0; // Force Calculation parameter 2

extern CProxy_Scheduler scheduler;

class GPUWork : public CBase_GPUWork {
public:
  GPUWork() {}
  GPUWork(CkMigrateMessage* msg) {}

  // L same row, U same col
  void gpu_offload(list<JMessage> msgs) {    
    //ckout << "---- offloading to GPU" << endl;

    char note[200];
    sprintf(note, "gpu offloading happening");
    traceUserSuppliedNote(note);

    vector<double> Lvec, Uvec, Avec;
    vector<int> Lstart, Lend, Ustart, Uend;

    int lastX = -1, lastY = -1;
    int startL = 0, startU = 0;

    int block = 0;

    for (list<JMessage>::iterator iter = msgs.begin(); 
	 iter != msgs.end(); ++iter) {
      bool copyL = true, copyU = true;
      int tsize = iter->fsize * iter->ssize;

      block = iter->fsize;

      if (lastX == iter->x)  
        copyL = false;
      else if (lastX != -1)
        startL += tsize;

      if (lastY == iter->y)
        copyU = false;
      else if (lastY != -1)
        startU += tsize;

      for (int i = 0; i < tsize; i++) {
        if (copyL)
          Lvec.push_back(iter->first[i]);

        if (copyU)
          Uvec.push_back(iter->second[i]);

        Avec.push_back(iter->LU[i]);

        Lstart.push_back(startL);
        Lend.push_back(startL + tsize);

        Ustart.push_back(startU);
        Uend.push_back(startU + tsize);
      }

      lastX = iter->x;
      lastY = iter->y;
    }

    double *Lm, *Um, *Am;
    int *Ls, *Le, *Us, *Ue;
    int size;

    Lm = &Lvec[0];
    Um = &Uvec[0];
    Am = &Avec[0];

    Ls = &Lstart[0];
    Le = &Lend[0];
    Us = &Ustart[0];
    Ue = &Uend[0];

    size = Avec.size();

    FakeGPUDGEMM(size, msgs.size(), block, Um, Lm, Am,
                 Us, Ue, Ls, Le);

    /*GPUinteractNew(Fx, Fy, Sx, Sy, Ffx, Ffy, Sfx, Sfy, fn, sn, 
      fS, fE, sS, sE);*/

    int firstLoc = 0, secondLoc = 0;

    for (list<JMessage>::iterator iter = msgs.begin();
	 iter != msgs.end(); ++iter) {
      int tsize = iter->fsize * iter->ssize;

      for (int i = 0; i < tsize; i++) {
        iter->LU[i] = Am[i + firstLoc];
      }

      firstLoc += tsize;
    }

    scheduler[CkMyPe()].finishedGPU(msgs);
  }

  void FakeGPUDGEMM(int tsize, int agglom, int block, double *Ap, double *Bp,
                    double *Cp, int *Astart, int *Aend, int *Bstart,
                    int *Bend) {
    int bsize = tsize / agglom;

    for (int m = 0; m < agglom; m++) {
      /*ckout << "bstart = " << Bstart[bsize * m] << endl;
        ckout << "astart = " << Astart[bsize * m] << endl;*/

      double *A = Astart[bsize * m] + Ap;
      double *B = Bstart[bsize * m] + Bp;
      double *C = m * bsize + Cp;

      /*for (int i = 0; i < block; i++) {
        ckout << "L = " << B[i] << endl;
        ckout << "U = " << A[i] << endl;
        ckout << "LU = " << C[i] << endl;
        }*/

      int N = block;
      unsigned i, j, k;

      for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) {
          const double *Ai_ = A + i;
          const double *B_j = B + j*N;

          double cij = *(C + j*N + i);

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
