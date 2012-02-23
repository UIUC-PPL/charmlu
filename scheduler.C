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

BlockScheduler::BlockScheduler(CProxy_LUBlk luArr_, LUConfig config, CProxy_LUMgr mgr_,
                                CProxy_StaticBlockSchedule staticProxy_, int numBlks_)
  : luArr(luArr_)
  , mgr(mgr_.ckLocalBranch())
  , maxMemory(0)
  , maxMemoryIncreases(0)
  , maxMemoryStep(-1)
  , numBlks(numBlks_)
  , nextRelease(-1) {
  staticProxy = staticProxy_;
  contribute(CkCallback(CkIndex_LUBlk::schedulerReady(NULL), luArr));
}

void BlockScheduler::registerBlock(CkIndex2D index) {
  //CkPrintf("registerBlock (%d,%d)\n", index.x, index.y);
  myBlocks.insert(index.x*numBlks+index.y);
  if (nextRelease == index.x * numBlks + index.y) {
    luArr(index.x,index.y).ckLocal()->releaseDep(0);
    nextRelease = -1;
  }
}

void BlockScheduler::unRegisterBlock(CkIndex2D index) {
  //CkPrintf("unRegisterBlock (%d,%d)\n", index.x, index.y);
  //fflush(stdout);
  assert(myBlocks.find(index.x*numBlks+index.y) != myBlocks.end());
  myBlocks.erase(index.x*numBlks+index.y);
}

void BlockScheduler::release(int x, int y) {
  if (myBlocks.find(x * numBlks + y) != myBlocks.end()) {
    luArr(x,y).ckLocal()->releaseDep(0);
  } else {
    nextRelease = x * numBlks + y;
  }
}

void BlockScheduler::notifyMigrate(int x, int y, int nproc) {
  luArr(x,y).ckLocal()->doMigrate(nproc);
}

void BlockScheduler::storeMsg(blkMsg* m) {
  CkIndex2D indx = m->indx;
  int step = CkGetRefNum(m);
  //CkPrintf("storeMsg from = (%d,%d), step = %d\n", m->indx.x, m->indx.y, step);
  //fflush(stdout);
  if(m->rightward) {
    for (int index = indx.y + 1; index < numBlks; index++) {  
        int proc = staticProxy.ckLocalBranch()->blockToProcs[indx.x * numBlks + index][step];
        if(CkMyPe() == proc) {
          CmiReference(UsrToEnv(m));
          CkIndex2D sIndx;
          sIndx.x = indx.x;
          sIndx.y = index;
          tryDeliver(m, sIndx);
        }
    }
    CmiFree(UsrToEnv(m));
  } else {
    for (int index = indx.x + 1; index < numBlks; index++) {
        int proc = staticProxy.ckLocalBranch()->blockToProcs[index * numBlks + indx.y][step];
        if(CkMyPe() == proc) {
          CmiReference(UsrToEnv(m));
          CkIndex2D sIndx;
          sIndx.x = index;
          sIndx.y = indx.y;
          tryDeliver(m, sIndx);
        }
    }
    CmiFree(UsrToEnv(m));
  }
}

void BlockScheduler::tryDeliver(blkMsg* m, CkIndex2D indx) {
  if(myBlocks.find(indx.x*numBlks+indx.y) != myBlocks.end()) {
    deliver(m,indx);
  }
  else {
    msgs[indx.x*numBlks+indx.y].push_back(m);  
  }
}

void BlockScheduler::deliver(blkMsg* m, CkIndex2D indx) {
  if(m->rightward) {
    if(m->indx.x == m->indx.y)
      luArr[indx].ckLocal()->recvL(m);
    else
      luArr[indx].ckLocal()->recvTrailingL(m);
  } else {
    if(m->indx.x == m->indx.y)
      luArr[indx].ckLocal()->recvU(m);
    else
      luArr[indx].ckLocal()->recvTrailingU(m);
  }
  //CmiFree
}

void BlockScheduler::checkMsgs(CkIndex2D indx) {
  for(list< blkMsg* >::iterator iter = msgs[indx.x*numBlks+indx.y].begin();
        iter != msgs[indx.x*numBlks+indx.y].end(); ++iter) {
    deliver(*iter, indx);
  }
  msgs[indx.x*numBlks+indx.y].clear();
}

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
  contribute(CkCallback(CkIndex_LUBlk::regDone(NULL), luArr));
}
