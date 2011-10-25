/**
 * @file A Charm++ implementation of LU
 */

#include "luConfig.h"

/// The build system should define this macro to be the commit identifier
#ifndef LU_REVISION
    #define LU_REVISION Unknown
#endif

#include <string.h>
#include <iostream>
#include <sstream>
#include <map>
#include <algorithm>
using std::min;
#include <limits>
#include <cmath>

#include "platformBlas.h"

#if USE_MEMALIGN
#include <malloc.h>
#endif

#include "lu.decl.h"
#include <trace-projections.h>
#include <ckmulticast.h>
#include <queueing.h> // for access to memory threshold setting

#include "lu.h"
#include "manager.h"
#include "mapping.h"
#include "messages.h"

// Define static variable
#if defined(LU_TRACING)
  int traceToggler::traceCmdHandlerID;
#endif

CkReductionMsg *maxMaxElm(int nMsg, CkReductionMsg **msgs) {
  CkAssert(nMsg > 0);

  MaxElm *l = (MaxElm*) msgs[0]->getData();
  for (int i = 1; i < nMsg; ++i) {
    MaxElm *n = (MaxElm *) msgs[i]->getData();
    if (fabs(n->val) > fabs(l->val))
      l = n;
  }

  return CkReductionMsg::buildNew(sizeof(MaxElm), l);
}

/// Global that holds the reducer type for MaxElm
CkReduction::reducerType MaxElmReducer;

/// Function that registers this reducer type on every processor
void registerMaxElmReducer() {
  MaxElmReducer = CkReduction::addReducer(maxMaxElm);
}

struct traceLU {
  int step, event;
  double startTime;
  traceLU(int internalStep, int eventType)
    : step(internalStep), event(eventType), startTime(CkWallTimer()) {
    traceUserSuppliedData(internalStep);
    traceMemoryUsage();
  }

  ~traceLU() {
    traceUserBracketEvent(event, startTime, CkWallTimer());
  }
};

void LUBlk::flushLogs() {
#if defined(LU_TRACING)
  flushTraceLog();
#endif
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
  suggestedPivotBatchSize = cfg.pivotBatchSize;
  schedAdaptMemThresholdMB = cfg.memThreshold;
  initDone = initialization;
  factorizationDone = factorization;
  solveDone = solution;
  internalStep = 0;

  CkAssert(blkSize > 0);

  traceUserSuppliedData(-1);
  traceMemoryUsage();

  CkMulticastMgr *mcastMgr = CProxy_CkMulticastMgr(cfg.mcastMgrGID).ckLocalBranch();

  /// Chares on the active panels will create sections of their brethren
#if defined(CHARMLU_USEG_FROM_BELOW)
  if (isOnDiagonal() || isBelowDiagonal())
#else
  if (isOnDiagonal())
#endif
  {
    // Elements in the active panel, not including this block
    CkVec<CkArrayIndex2D> activeElems;
    for (int i = thisIndex.y+1; i < numBlks; i++)
      if (i != thisIndex.x)
        activeElems.push_back(CkArrayIndex2D(i, thisIndex.y));
    activePanel = CProxySection_LUBlk::ckNew(thisArrayID, activeElems.getVec(), activeElems.size());
    activePanel.ckSectionDelegate(mcastMgr);
    rednSetupMsg *activePanelMsg = new rednSetupMsg(cfg.mcastMgrGID);
    activePanel.prepareForActivePanel(activePanelMsg);
  }
  /// Chares on the array diagonal will now create pivot sections that they will talk to
  if (isOnDiagonal()) {
    // Create the pivot section
    pivotSection = CProxySection_LUBlk::ckNew(thisArrayID, thisIndex.x,numBlks-1,1,thisIndex.y,thisIndex.y,1);
    rowBeforeDiag = CProxySection_LUBlk::ckNew(thisArrayID, thisIndex.x,thisIndex.x,1,0,thisIndex.y-1,1);
    rowAfterDiag = CProxySection_LUBlk::ckNew(thisArrayID, thisIndex.x,thisIndex.x,1,thisIndex.y+1,numBlks-1,1);
    // Delegate pivot section to the manager
    pivotSection.ckSectionDelegate(mcastMgr);
    rowBeforeDiag.ckSectionDelegate(mcastMgr);
    rowAfterDiag.ckSectionDelegate(mcastMgr);
    // Set the reduction client for this pivot section
    mcastMgr->setReductionClient( pivotSection, new CkCallback( CkReductionTarget(LUBlk,colMax), thisProxy(thisIndex.y, thisIndex.y) ) );

    // Invoke a dummy mcast so that all the section members know which section to reduce along
    rednSetupMsg *pivotMsg = new rednSetupMsg(cfg.mcastMgrGID);
    rednSetupMsg *rowBeforeMsg = new rednSetupMsg(cfg.mcastMgrGID);
    rednSetupMsg *rowAfterMsg = new rednSetupMsg(cfg.mcastMgrGID);

    pivotSection.prepareForPivotRedn(pivotMsg);
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
  delete LUmsg;
  LU = NULL;
}

void LUBlk::computeU(double *LMsg) {
  traceLU t(internalStep, cfg.traceComputeU);
#if USE_ESSL || USE_ACML
  // LMsg is implicitly transposed by telling dtrsm that it is a
  // right, upper matrix. Since this also switches the order of
  // multiplication, the transpose is output to LU.
  dtrsm(BLAS_RIGHT, BLAS_UPPER, BLAS_NOTRANSPOSE, BLAS_UNIT, blkSize, blkSize, 1.0, LMsg, blkSize, LU, blkSize);
#else
  cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, blkSize, blkSize, 1.0, LMsg, blkSize, LU, blkSize);
#endif
}

void LUBlk::updateMatrix(double *incomingL, double *incomingU) {
  traceLU t(internalStep, cfg.traceTrailingUpdate);

#if USE_ESSL || USE_ACML
  // By switching the order of incomingU and incomingL the transpose
  // is applied implicitly: C' = B*A
  dgemm( BLAS_NOTRANSPOSE, BLAS_NOTRANSPOSE,
         blkSize, blkSize, blkSize,
         -1.0, incomingU,
         blkSize, incomingL, blkSize,
         1.0, LU, blkSize);
#else
  cblas_dgemm( CblasRowMajor,
               CblasNoTrans, CblasNoTrans,
               blkSize, blkSize, blkSize,
               -1.0, incomingL,
               blkSize, incomingU, blkSize,
               1.0, LU, blkSize);
#endif
}

void LUBlk::setupMsg(bool reverse) {
  // Setup multicast of message to a dynamic set of processors
  blkMsg *m = LUmsg;

  CkAssert(requestingPEs.size() <= maxRequestingPEs);

  std::sort(requestingPEs.begin(), requestingPEs.end());
  if (reverse) std::reverse(requestingPEs.begin(), requestingPEs.end());

  DEBUG_PRINT("Preparing block for delivery to %d PEs", requestingPEs.size());

  // Junk value to catch bugs
  m->npes_sender = -1;
  m->npes_receiver = requestingPEs.size();
  m->offset = 0;
  memcpy(m->pes, &requestingPEs[0], sizeof(requestingPEs[0])*m->npes_receiver);

  requestingPEs.clear();
}

// Schedule U to be sent downward to the blocks in the same column
inline void LUBlk::scheduleDownwardU() {
  traceUserSuppliedData(internalStep);
  traceMemoryUsage();
  mgr->setPrio(LUmsg, MULT_RECV_U);

  DEBUG_PRINT("Multicast to part of column %d", thisIndex.y);
  localScheduler->scheduleSend(thisIndex, internalStep == thisIndex.y - 1);
}

// Schedule L to be sent rightward to the blocks in the same row
inline void LUBlk::scheduleRightwardL() {
  traceUserSuppliedData(internalStep);
  traceMemoryUsage();
  mgr->setPrio(LUmsg, MULT_RECV_L);

  DEBUG_PRINT("Multicast block to part of row %d", thisIndex.x);
  localScheduler->scheduleSend(thisIndex, true);
}

void LUBlk::requestBlock(int pe, int rx, int ry) {
  requestingPEs.push_back(pe);
  if (factored) {
    bool onActive = false;
    if      (isBelowDiagonal() && internalStep == thisIndex.y)
      onActive = true;
    else if (isAboveDiagonal() && internalStep == thisIndex.y - 1)
      onActive = true;

    localScheduler->scheduleSend(thisIndex, onActive);
  } else {
    DEBUG_PRINT("Queueing remote block for pe %d", pe);
  }
}

double* LUBlk::accessLocalBlock() {
  return LU;
}

void LUBlk::offDiagSolve(BVecMsg *m) {
  if (isOnDiagonal())
    return;

  // Do local portion of solve (daxpy)
  double *xvec = new double[blkSize], *preVec = m->data;
  for (int i = 0; i < blkSize; i++) {
    xvec[i] = 0.0;
    for (int j = 0; j < blkSize; j++)
      xvec[i] += LU[getIndex(i,j)] * preVec[j];
  }

  // Set the diagonal chare on my row as target of reduction
  CkCallback cb(CkIndex_LUBlk::recvSolveData(0), thisProxy(thisIndex.x, thisIndex.x));
  // Reduce row towards diagonal chare
  mcastMgr->contribute(sizeof(double) * blkSize, xvec, CkReduction::sum_double,
		       m->forward ? rowBeforeCookie : rowAfterCookie, cb, thisIndex.x);
  delete[] xvec;
}

// Copy received pivot data into its place in this block
void LUBlk::applySwap(int row, int offset, const double *data, double b) {
  DEBUG_PIVOT("(%d, %d): remote pivot inserted at %d\n", thisIndex.x, thisIndex.y, row);
  bvec[row] = b;
  memcpy( &(LU[getIndex(row,offset)]), data, sizeof(double)*(blkSize-offset) );
}

// Exchange local data
void LUBlk::swapLocal(int row1, int row2, int offset) {
  if (row1 == row2) return;
  std::swap(bvec[row1], bvec[row2]);
  /// @todo: Is this better or is it better to do 3 memcpys
  std::swap_ranges(&(LU[getIndex(row1,offset)]), &(LU[getIndex(row1,blkSize)]), &(LU[getIndex(row2,offset)]) );
}

void LUBlk::doPivotLocal(int row1, int row2) {
  // The chare indices where the two rows are located
  row1Index = row1 / blkSize;
  row2Index = row2 / blkSize;
  // The local indices of the two rows within their blocks
  localRow1 = row1 % blkSize;
  localRow2 = row2 % blkSize;
  remoteSwap = false;

  // If this block holds portions of both the current row and pivot row, its a local swap
  if (row1Index == thisIndex.x && row2Index == thisIndex.x) {
    swapLocal(localRow1, localRow2);
    // else if this block holds portions of at just one row, its a remote swap
  } else if (row1Index == thisIndex.x) {
    thisLocalRow = localRow1;
    otherRowIndex = row2Index;
    globalThisRow = row1;
    globalOtherRow = row2;
    remoteSwap = true;
  } else if (row2Index == thisIndex.x) {
    thisLocalRow = localRow2;
    otherRowIndex = row1Index;
    globalThisRow = row2;
    globalOtherRow = row1;
    remoteSwap = true;
  }
  // else this block has no data affected by this pivot op
}

/// Record the effect of a pivot operation in terms of actual row numbers
void LUBlk::recordPivot(const int r1, const int r2) {
  numRowsSinceLastPivotSend++;
  // If the two rows are the same, then dont record the pivot operation at all
  if (r1 == r2) return;
  std::map<int,int>::iterator itr1, itr2;
  // The records for the two rows (already existing or freshly created)
  itr1 = (pivotRecords.insert(std::make_pair(r1,r1))).first;
  itr2 = (pivotRecords.insert(std::make_pair(r2,r2))).first;
  // Swap the values (the actual rows living in these two positions)
  std::swap(itr1->second, itr2->second);
}

/**
 * Is it time to send out the next batch of pivots?
 * @note: Any runtime adaptivity should be plugged here
 */
bool LUBlk::shouldSendPivots() {
  return (numRowsSinceLastPivotSend >= suggestedPivotBatchSize);
}

/// Periodically send out the agglomerated pivot operations
void LUBlk::announceAgglomeratedPivots() {
#if defined(VERBOSE_PIVOT_RECORDING) || defined(VERBOSE_PIVOT_AGGLOM)
  std::stringstream pivotLog;
  pivotLog<<"["<<thisIndex.x<<","<<thisIndex.y<<"]"
          <<" announcing "<<pivotRecords.size()<<" pivot operations in batch "<<pivotBatchTag<<std::endl;
#endif

  // Create and initialize a msg to carry the pivot sequences
  pivotSequencesMsg *msg = new (numRowsSinceLastPivotSend+1, numRowsSinceLastPivotSend*2, sizeof(int)*8) pivotSequencesMsg(pivotBatchTag, numRowsSinceLastPivotSend);
  msg->numSequences = 0;
  memset(msg->seqIndex, 0, sizeof(int) * numRowsSinceLastPivotSend);
  memset(msg->pivotSequence, 0, sizeof(int) * numRowsSinceLastPivotSend*2);

  /// Parse the pivot operations and construct optimized pivot sequences
  int seqNo = -1, i = 0;
  std::map<int,int>::iterator itr = pivotRecords.begin();
  while (itr != pivotRecords.end()) {
#ifdef VERBOSE_PIVOT_RECORDING
    pivotLog<<std::endl;
#endif
    msg->seqIndex[++seqNo] = i;
    int chainStart = itr->first;
    msg->pivotSequence[i++] = chainStart;
#ifdef VERBOSE_PIVOT_RECORDING
    pivotLog<<chainStart;
#endif
    while (itr->second != chainStart) {
      msg->pivotSequence[i] = itr->second;
#ifdef VERBOSE_PIVOT_RECORDING
      pivotLog<<" <-- "<<itr->second;
#endif
      std::map<int,int>::iterator prev = itr;
      itr = pivotRecords.find(itr->second);
      pivotRecords.erase(prev);
      i++;
    }
    pivotRecords.erase(itr);
    itr = pivotRecords.begin();
  }
  msg->seqIndex[++seqNo] = i; ///< @note: Just so that we know where the last sequence ends
  msg->numSequences = seqNo;
#if defined(VERBOSE_PIVOT_RECORDING) || defined(VERBOSE_PIVOT_AGGLOM)
  CkPrintf("%s\n", pivotLog.str().c_str());
#endif

  mgr->setPrio(msg, PIVOT_RIGHT_SEC, -1, thisIndex.y);
  thisProxy.applyTrailingPivots(msg);

  // Prepare for the next batch of agglomeration
  pivotRecords.clear();
  pivotBatchTag += numRowsSinceLastPivotSend;
  numRowsSinceLastPivotSend = 0;
}

/// Given a set of pivot ops, send out participating row chunks that you own
void LUBlk::sendPendingPivots(const pivotSequencesMsg *msg) {
#ifdef VERBOSE_PIVOT_AGGLOM
  std::stringstream pivotLog;
  pivotLog << "[" << thisIndex.x << "," << thisIndex.y << "]"
           <<" processing "<< msg->numSequences
           <<" pivot sequences in batch " << pivotBatchTag;
#endif

  const int *pivotSequence = msg->pivotSequence, *idx = msg->seqIndex;
  int numSequences = msg->numSequences;

  // Count the number of rows that I send to each chare
  int numMsgsTo[numBlks];
  memset(numMsgsTo, 0, sizeof(int)*numBlks);
  for (int i = 0; i < numSequences; i++) {
    for (int j = idx[i]; j < idx[i+1]; j++) {
      int recverIdx = pivotSequence[j]/blkSize;
      int senderIdx = -1;
      // circular traversal of pivot sequence
      if (j < idx[i+1]-1)
        senderIdx = pivotSequence[j+1]/blkSize;
      else
        senderIdx = pivotSequence[idx[i]]/blkSize;
      // If this chare the sending to another chare
      if (thisIndex.x == senderIdx && thisIndex.x != recverIdx)
        numMsgsTo[recverIdx]++;
    }
  }

  // Preallocate messages that are large enough to carry the rows that this
  // block will send to each of the other chares
  pivotRowsMsg* outgoingPivotMsgs[numBlks];
  memset(outgoingPivotMsgs, 0, sizeof(pivotRowsMsg*) * numBlks);
  for (int i = 0; i < numBlks; i++) {
    // If this other chare expects row data from me
    if (numMsgsTo[i] > 0) {
      // Create a big enough msg to carry all the rows I'll be sending to this chare
      outgoingPivotMsgs[i] = new (numMsgsTo[i], numMsgsTo[i]*blkSize, numMsgsTo[i], sizeof(int)*8)
        pivotRowsMsg(blkSize, pivotBatchTag);
      // Set a priority thats a function of your location wrt to the critical path
      if (thisIndex.y < internalStep)
        mgr->setPrio(outgoingPivotMsgs[i], PIVOT_NOT_CRITICAL);
      else
        mgr->setPrio(outgoingPivotMsgs[i], PIVOT_CRITICAL);
    }
  }

  pendingIncomingPivots = 0;
  double *tmpBuf = 0, tmpB;

  // Parse each sequence independently
  for (int i=0; i < numSequences; i++) {
#ifdef VERBOSE_PIVOT_AGGLOM
    pivotLog << "\n[" << thisIndex.x << "," << thisIndex.y << "] sequence " << i << ": ";
#endif

    // Find the location of this sequence in the msg buffer
    const int *first = pivotSequence + idx[i];
    const int *beyondLast = pivotSequence + idx[i+1];
    CkAssert(beyondLast - first >= 2);

    // Identify a remote row in the pivot sequence as a point at which to
    // start and stop processing the circular pivot sequence
    const int *ringStart = first;
    while ((*ringStart / blkSize == thisIndex.x) && (ringStart < beyondLast))
      ringStart++;
    const int *ringStop = ringStart;

    // If there are no remote rows in the sequence, we *have* to use a tmp buffer
    // The tmp buffer will now complete the circular sequence
    bool isSequenceLocal = false;
    if (ringStart == beyondLast) {
      isSequenceLocal = true;
      ringStart = first;
      ringStop  = beyondLast - 1;

      if (tmpBuf == 0) tmpBuf = new double[blkSize];
      memcpy(tmpBuf, &LU[getIndex(*ringStart%blkSize,0)], blkSize*sizeof(double));
      tmpB = bvec[*ringStart%blkSize];
#ifdef VERBOSE_PIVOT_AGGLOM
      pivotLog << "tmp <-cpy- "<< *ringStart << "; ";
#endif
    }

    // Process all the pivot operations in the circular sequence
    const int *to = ringStart;
    do {
      const int *from = (to + 1 == beyondLast) ? first : to + 1;
      int fromChareIdx = *from / blkSize;
      // If the current source row in the pivot sequence belongs to me, send it
      if (fromChareIdx == thisIndex.x) {
        int toChareIdx = *to / blkSize;
        int fromLocal  = *from % blkSize;
        // If you're sending to yourself, memcopy
        if (toChareIdx == thisIndex.x) {
          applySwap(*to%blkSize, 0, &LU[getIndex(fromLocal,0)], bvec[fromLocal]);
#ifdef VERBOSE_PIVOT_AGGLOM
          pivotLog << *to << " <-cpy- " << *from << "; ";
#endif
        }
        // else, copy the data into the appropriate msg
        else {
          outgoingPivotMsgs[*to/blkSize]->copyRow(*to, &LU[getIndex(fromLocal,0)], bvec[fromLocal]);
#ifdef VERBOSE_PIVOT_AGGLOM
          pivotLog << *to << " <-msg- " << *from << "; ";
#endif
        }
      }
      // else, the source data is remote
      else {
        // if the current destination row belongs to me, make sure I expect the remote data
        if (*to / blkSize == thisIndex.x) {
          pendingIncomingPivots++;
#ifdef VERBOSE_PIVOT_AGGLOM
          pivotLog << *to << " <-inwd- " << *from << "; ";
#endif
        }
        // else, i dont worry about this portion of the exchange sequence which is completely remote
        else {
#ifdef VERBOSE_PIVOT_AGGLOM
          pivotLog << *to << " <-noop- " << *from << "; ";
#endif
        }
      }
      // Setup a circular traversal of the pivot sequence
      if (++to == beyondLast)
        to = first;
    } while (to != ringStop); // Keep going till you complete the ring

    // If the sequence was completely local, complete the circular sequence
    // by copying the temp buffer back into the matrix block
    if (isSequenceLocal) {
      applySwap(*(beyondLast - 1) % blkSize, 0, tmpBuf, tmpB);
#ifdef VERBOSE_PIVOT_AGGLOM
      pivotLog << *(beyondLast-1) << "<-cpy- tmp " << "; ";
#endif
    }
  } // end for loop through all sequences

  // Send out all the msgs carrying pivot data to other chares
  for (int i=0; i< numBlks; i++)
    if (numMsgsTo[i] > 0)
      thisProxy(i, thisIndex.y).trailingPivotRowsSwap(outgoingPivotMsgs[i]);

  if (tmpBuf) delete [] tmpBuf;
#ifdef VERBOSE_PIVOT_AGGLOM
  CkPrintf("%s\n",pivotLog.str().c_str());
#endif
}

// Internal functions for creating messages to encapsulate the priority
blkMsg* LUBlk::createABlkMsg() {
  int prioBits = mgr->bitsOfPrio();
  maxRequestingPEs = CProxy_LUMap(cfg.map).ckLocalBranch()->pesInPanel(thisIndex);
  blkMsg *msg = new (blkSize*blkSize, maxRequestingPEs, prioBits) blkMsg(thisIndex);
  memset(msg->pes, -1, maxRequestingPEs*sizeof(int));
  return msg;
}

MaxElm LUBlk::findMaxElm(int startRow, int col, MaxElm first) {
  MaxElm l = first;
  for (int row = startRow; row < blkSize; row++)
    if (fabs(LU[getIndex(row,col)]) > fabs(l.val)) {
      l.val = LU[getIndex(row,col)];
      l.loc = row + blkSize * thisIndex.x;
    }
  return l;
}

/// Update the sub-block of this L block starting at specified
/// offset from the active column
void LUBlk::updateLsubBlock(int activeCol, double* U, int offset, int startingRow) {
  // Should only get called on L blocks
  CkAssert(isOnDiagonal() || isBelowDiagonal());
  // Check for input edge cases
  if ((activeCol + offset) >= blkSize || startingRow >= blkSize)
    return;
#if USE_ESSL || USE_ACML
  dger(blkSize-(activeCol+offset), blkSize-startingRow,
       -1.0,
       U+offset, 1,
       &LU[getIndex(startingRow,activeCol)], blkSize,
       &LU[getIndex(startingRow,activeCol+offset)], blkSize);
#elif USE_ACCELERATE_BLAS
  for(int j = startingRow; j < blkSize; j++)
    for(int k = activeCol+offset; k<blkSize; k++)
      LU[getIndex(j,k)] -=  LU[getIndex(j,activeCol)] * U[k-activeCol];
#else
  cblas_dger(CblasRowMajor,
             blkSize-startingRow, blkSize-(activeCol+offset),
             -1.0,
             &LU[getIndex(startingRow,activeCol)], blkSize,
             U+offset, 1,
             &LU[getIndex(startingRow,activeCol+offset)], blkSize);
#endif
}


/// Compute the multipliers based on the pivot value in the
/// received row of U and also find the candidate pivot in
/// the immediate next column (after updating it simultaneously)
MaxElm LUBlk::computeMultipliersAndFindColMax(int col, double *U, int startingRow) {
  // Should only get called on L blocks
  CkAssert(isOnDiagonal() || isBelowDiagonal());
  MaxElm maxVal;
  // Check for input edge cases
  if (col >= blkSize || startingRow >= blkSize)
    return maxVal;

  if (col < blkSize -1) {
    for (int j = startingRow; j < blkSize; j++) {
      // Compute the multiplier
      LU[getIndex(j,col)]    = LU[getIndex(j,col)] / U[0];
      // Update the immediate next column
      LU[getIndex(j,col+1)] -= LU[getIndex(j,col)] * U[1];
      // Update the max value thus far
      if ( fabs(LU[getIndex(j,col+1)]) > fabs(maxVal.val) ) {
        maxVal.val = LU[getIndex(j,col+1)];
        maxVal.loc = j;
      }
    }
    // Convert local row num to global rownum
    maxVal.loc += thisIndex.x*blkSize;
  }
  else {
    for (int j = startingRow; j < blkSize; j++)
      LU[getIndex(j,col)]    = LU[getIndex(j,col)] / U[0];
  }

  return maxVal;
}

#include "luUtils.def.h"
#include "lu.def.h"
#include "luMessages.def.h"
