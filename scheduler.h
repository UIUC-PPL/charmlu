/**
 * BlockScheduler
 *
 * The block scheduler maintains a list of the blocks assigned to its
 * processor, and tracks what step they have reached. Within the bounds of the
 * memory threshold, it requests blocks from remote processors that are needed
 * for local triangular solves and trailing updates. To eliminate the
 * possibility of deadlock, the order in which operations are executed, and
 * hence remote blocks requested, must be carefully selected.
 *
 * Instead for forcing updates to execute in step-order, the order that blocks
 * can execute is selected based on their dependencies on other blocks. The
 * modeling of these dependencies and the demonstration that this leads to a
 * deadlock-free order is shown in a technical report, "Exploring Partial
 * Synchrony in an Asynchronous Environment Using Dense LU".
 */

#ifndef LU_SCHEDULER_H
#define LU_SCHEDULER_H

#include "messages.h"
#include "luConfig.h"
#include "lu.decl.h"
#include <list>
#include <map>
#include <set>
#include <utility>
#include <algorithm>
#include <utility>

/**
 * The BlockScheduler that tracks the dependencies of blocks, requests blocks,
 * and thereby manages memory consumption on a particular processor.
 */
class BlockScheduler : public CBase_BlockScheduler {
public:
  BlockScheduler(CProxy_LUBlk luArr_, LUConfig config, CProxy_LUMgr mgr_,
                  CProxy_StaticBlockSchedule, int);
  BlockScheduler(CkMigrateMessage *m) { }

  // All LUBlk chares that live on the same processor, must register with the
  // BlockScheduler at startup
  void registerBlock(CkIndex2D index);
  void unRegisterBlock(CkIndex2D index);
  // Reduction to indicate that all chares have registered
  void allRegistered(CkReductionMsg *m);
  void outputStats();
  void tryDeliver(blkMsg* m, CkIndex2D indx);
  void deliver(blkMsg* m, CkIndex2D indx);
  void storeMsg(blkMsg* m);
  void checkMsgs(CkIndex2D indx);
  void release(int x, int y);
  void notifyMigrate(int, int, int);
  void warmup(rednSetupMsg* msg);
  
private:
  LUMgr *mgr;
  CProxy_LUBlk luArr;
  set<int> myBlocks;
  CProxy_StaticBlockSchedule staticProxy;
  map< int, list< blkMsg* > > msgs;
  int numBlks;
  int nextRelease;

  // Memory usage instrumentation
  // Memory occupied before factorization starts (kilobytes)
  int baseMemory;
  // Maximum memory usage seen (kilobytes)
  int maxMemory, maxMemoryIncreases, maxMemoryStep;
};

#endif
