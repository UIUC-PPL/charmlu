/**
 * @file A Charm++ implementation of LU
 */

#include <algorithm>
using std::min;
#include <cmath>
#include <limits>
#include <map>
#include <string.h>

#include "lu.h"
#include "mapping.h"
#include "platformBlas.h"

// Execute a trailing update
void LUBlk::updateMatrix(double *incomingL, double *incomingU) {
#if USE_ESSL || USE_ACML
  // By switching the order of incomingU and incomingL the transpose
  // is applied implicitly: C' = B*A
  dgemm(BLAS_NOTRANSPOSE, BLAS_NOTRANSPOSE,
        blkSize, blkSize, blkSize,
        -1.0, incomingU,
        blkSize, incomingL, blkSize,
        1.0, LU, blkSize);
#else
  cblas_dgemm(CblasRowMajor,
              CblasNoTrans, CblasNoTrans,
              blkSize, blkSize, blkSize,
              -1.0, incomingL,
              blkSize, incomingU, blkSize,
              1.0, LU, blkSize);
#endif
}

void LUBlk::decompose() {
  for (int j = 0; j < blkSize; j++) {
    for (int i = 0; i <= j; i++) {
      double sum = 0.0;
      for (int k = 0; k < i; k++)
        sum += LU[getIndex(i, k)] * LU[getIndex(k, j)];
      LU[getIndex(i, j)] -= sum;
    }
    for (int i = j + 1; i < blkSize; i++) {
      double sum = 0.0;
      for (int k = 0; k < j; k++)
        sum += LU[getIndex(i, k)] * LU[getIndex(k, j)];
      LU[getIndex(i, j)] -= sum;
      LU[getIndex(i, j)] /= LU[getIndex(j, j)];
    }
  }
}

// Compute a triangular solve
void LUBlk::computeU(double *LMsg) {
#if USE_ESSL || USE_ACML
  // LMsg is implicitly transposed by telling dtrsm that it is a
  // right, upper matrix. Since this also switches the order of
  // multiplication, the transpose is output to LU.
  dtrsm(BLAS_RIGHT, BLAS_UPPER, BLAS_NOTRANSPOSE, BLAS_UNIT, blkSize, blkSize, 1.0, LMsg, blkSize, LU, blkSize);
#else
  cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, blkSize, blkSize, 1.0, LMsg, blkSize, LU, blkSize);
#endif
}

// Compute a triangular solve
void LUBlk::computeL(double *UMsg) {
#if USE_ESSL || USE_ACML
  dtrsm(BLAS_LEFT, BLAS_LOWER, BLAS_NOTRANSPOSE, BLAS_UNIT, blkSize, blkSize, 1.0, UMsg, blkSize, LU, blkSize);
#else
  cblas_dtrsm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasUnit, blkSize, blkSize, 1.0, UMsg, blkSize, LU, blkSize);
#endif
}

// Send U downward to the blocks in the same column
inline void LUBlk::sendDownwardU() {
  blkMsg* msg = createABlkMsg();
  //mgr->setPrio(msg, MULT_RECV_U);
  CProxySection_LUBlk col = CProxySection_LUBlk::ckNew(thisArrayID,
                                                       thisIndex.x + 1, numBlks - 1, 1,
                                                       thisIndex.y, thisIndex.y, 1);
  col.ckSectionDelegate(mcastMgr);
  if (isOnDiagonal()) {
    col.recvU(msg);
  } else {
    col.recvTrailingU(msg);
  }
}

// Send L rightward to the blocks in the same row
inline void LUBlk::sendRightwardL() {
  blkMsg* msg = createABlkMsg();
  //mgr->setPrio(msg, MULT_RECV_L);
  CProxySection_LUBlk row = CProxySection_LUBlk::ckNew(thisArrayID, thisIndex.x,
                                                       thisIndex.x, 1, thisIndex.y + 1,
                                                       numBlks - 1, 1);
  row.ckSectionDelegate(mcastMgr);
  if (isOnDiagonal()) {
    row.recvL(msg);
  } else {
    row.recvTrailingL(msg);
  }
}

// Internal functions for creating messages to encapsulate the priority
blkMsg* LUBlk::createABlkMsg() {
  int prioBits = mgr->bitsOfPrio();
  blkMsg* msg = new (blkSize * blkSize, prioBits)
    blkMsg(LU, std::min(thisIndex.x, thisIndex.y), blkSize);
  return msg;
}

void LUBlk::init(const LUConfig _cfg, CProxy_LUMgr _mgr,
                 CProxy_BlockScheduler bs,
		 CkCallback initialization, CkCallback factorization, CkCallback solution) {
  scheduler = bs;
  localScheduler = scheduler[CkMyPe()].ckLocal();
  CkAssert(localScheduler);
  localScheduler->registerBlock(thisIndex);
  contribute(CkCallback(CkIndex_BlockScheduler::allRegistered(NULL), bs));
  cfg = _cfg;
  blkSize = cfg.blockSize;
  numBlks = cfg.numBlocks;
  mgr = _mgr.ckLocalBranch();
  initDone = initialization;
  factorizationDone = factorization;
  solveDone = solution;
  internalStep = 0;

  CkAssert(blkSize > 0);

  mcastMgr = CProxy_CkMulticastMgr(cfg.mcastMgrGID).ckLocalBranch();

  /// Chares on the active panels will create sections of their brethren
  if (isOnDiagonal()) {
    activePanel = CProxySection_LUBlk::ckNew(thisArrayID,
                                             thisIndex.x + 1, numBlks - 1, 1,
                                             thisIndex.y, thisIndex.y, 1);
    activePanel.ckSectionDelegate(mcastMgr);
    rednSetupMsg *activePanelMsg = new rednSetupMsg(cfg.mcastMgrGID);
    activePanel.prepareForActivePanel(activePanelMsg);
  }
  /// Chares on the array diagonal will now create pivot sections that they will talk to
  if (isOnDiagonal()) {
    rowBeforeDiag = CProxySection_LUBlk::ckNew(thisArrayID,
                                               thisIndex.x, thisIndex.x, 1, 0,
                                               thisIndex.y - 1, 1);
    rowAfterDiag = CProxySection_LUBlk::ckNew(thisArrayID, thisIndex.x,
                                              thisIndex.x, 1, thisIndex.y + 1,
                                              numBlks - 1, 1);
    rowBeforeDiag.ckSectionDelegate(mcastMgr);
    rowAfterDiag.ckSectionDelegate(mcastMgr);

    // Invoke a dummy mcast so that all the section members know which section to reduce along
    rednSetupMsg *rowBeforeMsg = new rednSetupMsg(cfg.mcastMgrGID);
    rednSetupMsg *rowAfterMsg = new rednSetupMsg(cfg.mcastMgrGID);

    rowBeforeDiag.prepareForRowBeforeDiag(rowBeforeMsg);
    rowAfterDiag.prepareForRowAfterDiag(rowAfterMsg);

    if (thisIndex.x == 0) {
      thisProxy.multicastRedns(0);
    }
  }
  // All chares except members of pivot sections are done with init
}

void LUBlk::prepareForActivePanel(rednSetupMsg *msg) { }

LUBlk::~LUBlk() {
  delete LU;
  LU = NULL;
}

#include "luUtils.def.h"
#include "lu.def.h"
#include "luMessages.def.h"
