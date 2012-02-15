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
using std::min;

/**
 * 2D chare array that embodies a block of the input matrix
 *
 * Main participant in the LU factorization and solve phases
 */
class LUBlk: public CBase_LUBlk {
public:
  /// Constructor
  LUBlk() {
    // allow SDAG to initialize its internal state for this chare
    __sdag_init();
  }

  void init(const LUConfig _cfg, CProxy_LUMgr _mgr, CProxy_BlockScheduler bs,
	    CkCallback initDone, CkCallback fznDone, CkCallback slnDone);
  void prepareForActivePanel(rednSetupMsg *msg);
  ~LUBlk();
  LUBlk(CkMigrateMessage* m) { CkAbort("LU blocks not migratable yet"); }
  // Added for migration
  void pup(PUP::er &p) {
    __sdag_pup(p);

    p | scheduler;
    p | cfg;
    p | blkSize;
    p | numBlks;
    p | startTime;
    p | activePanel;
    p | rowBeforeDiag;
    p | rowAfterDiag;
    p | rowBeforeCookie;
    p | rowAfterCookie;
    p | internalStep;
    p | updateExecuted;
    p | initDone; 
    p | factorizationDone;
    p | solveDone;
    p | mgrp;

    PUParray(p, LU, blkSize * blkSize);

    if (p.isUnpacking()) {
      mcastMgr = CProxy_CkMulticastMgr(cfg.mcastMgrGID).ckLocalBranch();
      mgr = mgrp.ckLocalBranch();
    }
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
