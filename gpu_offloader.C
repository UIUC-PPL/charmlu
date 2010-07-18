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
    ckout << "---- offloading to GPU" << endl;

    char note[200];
    sprintf(note, "gpu offloading happening");
    traceUserSuppliedNote(note);

    vector<double> Lvec, Uvec, Avec;
    vector<int> Lstart, Lend, Ustart, Uend;

    int lastX = -1, lastY = -1;
    int startL = 0, startU = 0;

    for (list<JMessage>::iterator iter = msgs.begin(); 
	 iter != msgs.end(); ++iter) {
      bool copyL = true, copyU = true;
      int tsize = iter->fsize * iter->ssize;

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

    /*GPUinteractNew(Fx, Fy, Sx, Sy, Ffx, Ffy, Sfx, Sfy, fn, sn, 
      fS, fE, sS, sE);*/

    /*firstLoc = 0;
      secondLoc = 0;*/

    /*for (list<JMessage>::iterator iter = msgs.begin(); 
	 iter != msgs.end(); ++iter) {
      
      for (int i = 0; i < iter->first->size(); i++) {
	(*(iter->first))[i].fx = Ffx[firstLoc];
	(*(iter->first))[i].fy = Ffy[firstLoc];
	firstLoc++;
      }

      for (int i = 0; i < iter->second->size(); i++) {
	(*(iter->second))[i].fx = Sfx[secondLoc];
	(*(iter->second))[i].fy = Sfy[secondLoc];
	secondLoc++;
      }
      }*/

    scheduler[CkMyPe()].finishedGPU(msgs);
  }

};


#include "gpuwork.def.h"

