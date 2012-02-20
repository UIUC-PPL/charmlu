#include "manager.h"
#include "messages.h"
#include "scheduler.h"
#include "luConfig.h"
#include "lu.decl.h"
#include <ckmulticast.h>
#include <limits>
#include <map>
#include <algorithm>
#include <vector>
#include <cassert>
using std::min;

struct StateNode {
  int taskid, taskType, is, proc;
};

struct StaticBlockSchedule : public CBase_StaticBlockSchedule {
  std::vector<std::vector<StateNode> > blockStates;
  int numBlks;
  
  StaticBlockSchedule(int numBlks)
    : numBlks(numBlks) {
    blockStates.resize(numBlks * numBlks);

    if (CkMyPe() == 0) {
      CkPrintf("reading schedule\n");
      fflush(stdout);
    }

    // read from file;
    int curBlock = -1, curVal = 0, x;
    FILE* schedule = fopen("schedule.dag", "r");
    StateNode node;
    while (fscanf(schedule, "%d ", &x) == 1) {
      curVal++;
      if (x == -1) {
	//CkPrintf("x == -1\n");
	curBlock++;
	curVal = 0;
      } else {
	switch (curVal) {
	case 1: node.taskid = x; break;
	case 2: node.taskType = x; break;
	case 3: node.is = x; break;
	case 4:
	  node.proc = x;
	  assert(curBlock < blockStates.size());
	  //CkPrintf("readNode %d: %d %d %d %d\n",
	  //curBlock, node.taskid, node.taskType, node.is, node.proc);
	  fflush(stdout);
	  blockStates[curBlock].push_back(node);
	  curVal = 0;
	  break;
	}
      }
    }
    fclose(schedule);

    if (CkMyPe() == 0) {
      CkPrintf("finished reading schedule\n");
      fflush(stdout);
    }
  }

  StateNode getNextState(int x, int y, int state) {
    int index = x * numBlks + y;
    assert(index < blockStates.size());
    
    //CkPrintf("(%d,%d): getNextState on index %d, current = %d\n", x, y, index, state);
    fflush(stdout);

    if (state + 1 == blockStates[index].size()) {
      StateNode state;
      state.taskType = -1;
      state.taskid = -1;
      state.is = -1;
      state.proc = -1;
      return state;
    } else {
      return blockStates[index][state + 1];
    }
  }
};

/**
 * 2D chare array that embodies a block of the input matrix
 *
 * Main participant in the LU factorization and solve phases
 */
class LUBlk: public CBase_LUBlk {
public:
  /// Constructor
  LUBlk() 
  : started(false)
  , currentState(-1) {
    // allow SDAG to initialize its internal state for this chare
    __sdag_init();
  }

  void init(const LUConfig _cfg, CProxy_LUMgr _mgr, CProxy_BlockScheduler bs,
	    CkCallback initDone, CkCallback fznDone, CkCallback slnDone,
	    CProxy_StaticBlockSchedule staticProxy_);
  void prepareForActivePanel(rednSetupMsg *msg);
  ~LUBlk();
  LUBlk(CkMigrateMessage* m) { CkAbort("LU blocks not migratable yet"); }
  // Added for migration
  void pup(PUP::er &p) {
    //CkPrintf("puping object unpacking = %s\n", p.isUnpacking() ? "true" : "false");

    p | blkSize;

    if (p.isUnpacking()) {
      LU = new double[blkSize * blkSize];
      bvec = new double[blkSize];
    }

    PUParray(p, LU, blkSize * blkSize);

    p | nproc;
    p | started;
    p | staticProxy;
    p | scheduler;
    p | cfg;
    p | numBlks;
    p | startTime;
    p | activePanel;
    p | rowBeforeDiag;
    p | rowAfterDiag;
    p | rowBeforeCookie;
    p | rowAfterCookie;
    p | currentState;
    p | internalStep;
    p | updateExecuted;
    p | initDone; 
    p | factorizationDone;
    p | solveDone;
    p | mgrp;

    if (p.isUnpacking()) {
      mcastMgr = CProxy_CkMulticastMgr(cfg.mcastMgrGID).ckLocalBranch();
      mgr = mgrp.ckLocalBranch();
    }

    __sdag_pup(p);

    if (p.isUnpacking()) {
      //CkPrintf("(%d,%d): migrated: internalStep = %d\n", thisIndex.x, thisIndex.y, internalStep);
      //CkPrintf("(%d,%d): migrated: numBlks = %d\n", thisIndex.x, thisIndex.y, numBlks);
      //CkPrintf("(%d,%d): migrated: blkSize = %d\n", thisIndex.x, thisIndex.y, blkSize);
    }
    fflush(stdout);
  }

  virtual void ckJustMigrated() {
    ArrayElement::ckJustMigrated();
    //CkPrintf("ckJustMigrated()\n");
    fflush(stdout);
    thisProxy[thisIndex].migrateDone(0);
  }

  void doMigrate(int proc) {
    migrateMe(proc);
  }

public:
  int internalStep;

protected:
  // Internal functions for creating messages to encapsulate the priority
  blkMsg* createABlkMsg();
  CProxy_BlockScheduler scheduler;

  /// Configuration settings
  LUConfig cfg;

  /// Variables used during factorization
  double *LU;
  int blkSize, numBlks;
  blkMsg *U, *L;

  int nproc;

  bool started;

  CProxy_StaticBlockSchedule staticProxy;
  CProxy_LUMgr mgrp;
  LUMgr *mgr;

  /// Variables used only during solution
  double *bvec;

  // Timer for each block
  double startTime;

  // The panel of blocks below the active diagonal chare
  CProxySection_LUBlk activePanel;
  CProxySection_LUBlk rowBeforeDiag;
  CProxySection_LUBlk rowAfterDiag;
  CkSectionInfo rowBeforeCookie;
  CkSectionInfo rowAfterCookie;

  int currentState;

  /// A pointer to the local branch of the multicast manager group that handles the pivot section comm
  CkMulticastMgr *mcastMgr;

  bool updateExecuted;
  CkCallback initDone, factorizationDone, solveDone;

  // System-generated macro for SDAG code
  LUBlk_SDAG_CODE

private:
  void decompose();
  /// Perform trailing update based on input matrices (dgemm)
  void updateMatrix(double *incomingL, double *incomingU);
  /// Performs the triangular solve required to compute a block of the U matrix (dtrsm)
  void computeU(double *LMsg);
  /// Performs the triangular solve required to compute a block of the L matrix (dtrsm)
  void computeL(double *UMsg);
  /// Send U downward to the blocks in the same column
  void sendDownwardU();
  // Send L rightwards to the blocks in the same row
  void sendRightwardL();

  inline int getIndex(int i, int j) {
    return i * blkSize + j;
  }

  bool isOnDiagonal()    { return thisIndex.x == thisIndex.y; }
  bool isAboveDiagonal() { return thisIndex.x <  thisIndex.y; }
  bool isBelowDiagonal() { return thisIndex.x >  thisIndex.y; }
};
