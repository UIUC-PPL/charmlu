/**
 * BlockScheduler
 *
 * Plan on a set of trailing updates up to the block (U or L) limit, unless the
 * block is local or shared by another trailing update. A block is eligible to
 * be planned if its necessary dependencies on other blocks have been
 * fulfilled, by that block having been previously planned. To simplify the
 * dependencies the block are split into regions (called panels in the code) to
 * signify that they have the same number of pending dependencies because they
 * are in the same region. This is explained further in a technical report,
 * "Exploring Partial Synchrony in an Asynchronous Environment Using Dense LU".
 */

#include "scheduler.h"
#include "lu.decl.h"
#include "messages.h"
#include "lu.h"
#include "register.h"

BlockScheduler::BlockScheduler(CProxy_LUBlk luArr_, LUConfig config, CProxy_LUMgr mgr_)
  : luArr(luArr_)
  , mgr(mgr_.ckLocalBranch())
  , maxMemory(0)
  , maxMemoryIncreases(0)
  , maxMemoryStep(-1) {
  contribute(CkCallback(CkIndex_LUBlk::schedulerReady(NULL), luArr));
}

void BlockScheduler::registerBlock(CkIndex2D index) { }

void printMemory(void *time, void *msg) {
  int *s = (int*) ((CkReductionMsg *)msg)->getData();
  int mem = s[0];
  CkPrintf("%s memory usage: %zu KiB, additional stats: ", time, mem);
  for (int i = 1; i < ((CkReductionMsg*)msg)->getLength() / sizeof(int); i++)
    CkPrintf("%zd ", s[i]);
  CkPrintf("\n");
  delete (CkReductionMsg*)msg;
}

void BlockScheduler::outputStats() {
  int stats[3];
  stats[0] = maxMemory;
  stats[1] = maxMemoryIncreases;
  stats[2] = maxMemoryStep;
  contribute(3 * sizeof(int), stats, CkReduction::max_int,
	     CkCallback(&printMemory, const_cast<char*>("Peak")));
}

void BlockScheduler::allRegistered(CkReductionMsg *m) {
  delete m;

  baseMemory = CmiMemoryUsage()/1024;
  contribute(sizeof(int), &baseMemory, CkReduction::max_int,
	     CkCallback(&printMemory, const_cast<char*>("Base")));
}
