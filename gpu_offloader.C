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

  void gpu_offload(list<JMessage> msgs) {    
    ckout << "---- offloading to GPU" << endl;

    char note[200];
    sprintf(note, "gpu offloading happening");
    traceUserSuppliedNote(note);

    float *Fx, *Fy, *Sx, *Sy, *Ffx, *Ffy, *Sfx, *Sfy;
    int *fS, *fE, *sS, *sE;
    int fn, sn;

    bool first = true;

    int firstLoc = 0, secondLoc = 0;

    vector<float> firstVecX, firstVecY, firstVecFX, firstVecFY;
    vector<float> secondVecX, secondVecY, secondVecFX, secondVecFY;
    vector<int> fStart, fEnd, sStart, sEnd;
    vector<int> fStarts, fEnds, sStarts, sEnds;

    for (list<JMessage>::iterator iter = msgs.begin(); 
	 iter != msgs.end(); ++iter) {
      
      int fstartLoc = firstLoc;
      int fstopLoc = firstLoc + iter->fsize - 1;
      int sstartLoc = secondLoc;
      int sstopLoc = secondLoc + iter->ssize - 1;

      fStarts.push_back(fstartLoc);
      fEnds.push_back(fstopLoc);
      sStarts.push_back(sstartLoc);
      sEnds.push_back(sstopLoc);
    }

    int ii = 0;

    for (list<JMessage>::iterator iter = msgs.begin(); 
	 iter != msgs.end(); ++iter) {
      
      int fstartLoc = firstLoc;
      int fstopLoc = firstLoc + iter->fsize - 1;
      int sstartLoc = secondLoc;
      int sstopLoc = secondLoc + iter->ssize - 1;

      for (int i = 0; i < iter->fsize; i++) {
	/*firstVecX.push_back((*(iter->first))[i].x);
	firstVecY.push_back((*(iter->first))[i].y);
	firstVecFX.push_back((*(iter->first))[i].fx);
	firstVecFY.push_back((*(iter->first))[i].fy);
	
	fStart.push_back(sStarts[ii]);
	fEnd.push_back(sEnds[ii]);*/

	/*assert(fstartLoc <= fstopLoc );
	assert(fstartLoc < MAX_OFFLOAD_PARTICLES);
	assert(fstopLoc < MAX_OFFLOAD_PARTICLES);*/

	firstLoc++;
      }

      for (int i = 0; i < iter->ssize; i++) {
	/*secondVecX.push_back((*(iter->second))[i].x);
	secondVecY.push_back((*(iter->second))[i].y);
	secondVecFX.push_back((*(iter->second))[i].fx);
	secondVecFY.push_back((*(iter->second))[i].fy);
	
	sStart.push_back(fStarts[ii]);
	sEnd.push_back(fEnds[ii]);*/

	/*assert(sstartLoc <= sstopLoc );
	assert(sstartLoc < MAX_OFFLOAD_PARTICLES);
	assert(sstopLoc < MAX_OFFLOAD_PARTICLES);*/

	secondLoc++;
      }

      ii++;
    }

    Fx = &firstVecX[0];
    Fy = &firstVecY[0];
    Sx = &secondVecX[0];
    Sy = &secondVecY[0];

    Ffx = &firstVecFX[0];
    Ffy = &firstVecFY[0];
    Sfx = &secondVecFX[0];
    Sfy = &secondVecFY[0];

    fS = &fStart[0];
    fE = &fEnd[0];
    sS = &sStart[0];
    sE = &sEnd[0];

    fn = firstVecX.size();
    sn = secondVecX.size();

    /*assert(fn == firstVecFX.size());
    assert(fn == firstVecFY.size());
    assert(sn == secondVecFX.size());
    assert(sn == secondVecFY.size());
    assert(fn == fStart.size());
    assert(fn == fEnd.size());
    assert(sn == sStart.size());
    assert(sn == sEnd.size());

    for (int i = 0; i < fn; i++) {
      assert(fStart[i] >= 0);
      assert(fStart[i] < fEnd[i]);
      assert(fEnd[i] < sn);
    }

    for (int i = 0; i < sn; i++) {
      assert(sStart[i] >= 0);
      assert(sStart[i] < sEnd[i]);
      assert(sEnd[i] < fn);
      }*/

    /*GPUinteractNew(Fx, Fy, Sx, Sy, Ffx, Ffy, Sfx, Sfy, fn, sn, 
      fS, fE, sS, sE);*/

    firstLoc = 0;
    secondLoc = 0;

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

