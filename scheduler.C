#include "scheduler.h"
#include "lu.decl.h"
#include "messages.h"
#include "lu.h"
#include <algorithm>

inline bool operator==(const CkIndex2D &l, const CkIndex2D &r)
{ return l.x == r.x && l.y == r.y; }
inline bool operator<(const CkIndex2D &l, const CkIndex2D &r)
{ return l.x < r.x || (l.x == r.x && l.y < r.y); }

void BlockScheduler::registerBlock(CkIndex2D index) {
  blockLimit--;
  if (index.x != 0 && index.y != 0)
    localBlocks.push_back(BlockState(index));
  CkAssert(blockLimit >= 2);
}

struct findByIndex {
  bool operator()(const BlockState &state) {
    return state.ix == index.x && state.iy == index.y;
  }

  CkIndex2D index;
  findByIndex(CkIndex2D _index) : index(_index) {}
};

StateList::iterator BlockScheduler::findBlockState(CkIndex2D index) {
  StateList::iterator iter = find_if(localBlocks.begin(),
                                     localBlocks.end(), findByIndex(index));
  if (iter == localBlocks.end()) {
    iter = find_if(pendingBlocks.begin(), pendingBlocks.end(), findByIndex(index));
  }
  if (iter == pendingBlocks.end()) {
    iter = find_if(readyBlocks.begin(), readyBlocks.end(), findByIndex(index));
  }
  CkAssert(iter != readyBlocks.end());
  return iter;
}

void BlockScheduler::pivotsDone(CkIndex2D index) {
  findBlockState(index)->pivotsDone = true;
  progress();
}

void BlockScheduler::dataReady(CkIndex2D index, BlockReadyMsg *m) {
  StateList::iterator iter = findBlockState(index);
  BlockState::InputState &input = m->src.x == index.x ? iter->Lstate : iter->Ustate;
  input.m = m;
  DEBUG_SCHED("dataReady message src = (%d, %d)\n", input.m->src.x, input.m->src.y);
  progress();
}

void BlockScheduler::updateDone(CkIndex2D index) {
  StateList::iterator iter = findBlockState(index);
  iter->updatesCompleted++;
  iter->pivotsDone = false;
  iter->Lstate.reset();
  iter->Ustate.reset();
}

struct updatesCompletedSorter {
  bool operator()(const BlockState& state1, const BlockState& state2) {
    return state1.updatesCompleted < state2.updatesCompleted;
  }
};

void BlockScheduler::progress() {
  bool stateModified = false;
  localBlocks.sort(updatesCompletedSorter());

  while ((pendingBlocks.size() + readyBlocks.size()) * 2 < blockLimit && localBlocks.size() > 0) {
    // Move some blocks from localBlocks to pendingBlocks
    pendingBlocks.splice(pendingBlocks.end(), localBlocks, localBlocks.begin());
    stateModified = true;
  }

  // Try to advance state for each
  for (StateList::iterator iter = pendingBlocks.begin(); iter != pendingBlocks.end();
       ++iter) {
    BlockState &block = *iter;
    DEBUG_SCHED("examining pending block (%d, %d)\n", block.ix, block.iy);
    if (block.pivotsDone) {
      // For Lstate
      if (block.Lstate.state == PENDING_SPACE) {
        block.Lstate.state = ALLOCATED;
        stateModified = true;
      }

      if (block.Lstate.state == ALLOCATED) {
        if (block.Lstate.m) {
          getBlock(block.Lstate);
          stateModified = true;
        } // else we wait for ready msg
      }

      // For Ustate
      if (block.Ustate.state == PENDING_SPACE) {
        block.Ustate.state = ALLOCATED;
        stateModified = true;
      }

      if (block.Ustate.state == ALLOCATED) {
        if (block.Ustate.m) {
          getBlock(block.Ustate);
          stateModified = true;
        } // else we wait for ready msg
      }

      if (block.Ustate.state == ARRIVED &&
          block.Lstate.state == ARRIVED) {
        DEBUG_SCHED("both block indicate arrived for (%d, %d)\n", block.ix, block.iy);
        readyBlocks.push_back(block);
        iter = pendingBlocks.erase(iter);
        stateModified = true;
      }
    }
  }

  StateList::iterator block = readyBlocks.begin();
  if (block != readyBlocks.end()) {
    CkIndex2D Lsrc = block->Lstate.m->src, Usrc = block->Ustate.m->src;
    double *Ldata = block->Lstate.data, *Udata = block->Ustate.data;
    CkAssert(block->updatesCompleted == Lsrc.y && block->updatesCompleted == Usrc.x);
    luArr(block->ix, block->iy).ckLocal()->
      processTrailingUpdate(block->updatesCompleted, (intptr_t)block->Lstate.data,
                            (intptr_t)block->Ustate.data);
    for (std::list<blkMsg*>::iterator iter = blockMessages.begin();
         iter != blockMessages.end();
         ++iter) {
      if ((*iter)->data == Ldata || (*iter)->data == Udata) {
        delete *iter;
        iter = blockMessages.erase(iter);
      }
    }
    if (block->updatesCompleted < std::min(block->ix, block->iy))
      localBlocks.splice(localBlocks.end(), readyBlocks, block);
    else
      readyBlocks.erase(block);
    stateModified = true;
  }

  if (stateModified) {
    progress();
  }
}

void BlockScheduler::getBlock(BlockState::InputState &input) {
  CkAssert(input.state == ALLOCATED);
  DEBUG_SCHED("requesting getBlock from (%d, %d)\n", input.m->src.x, input.m->src.y);
  luArr(input.m->src.x, input.m->src.y).
    getBlock(CkCallback(CkIndex_BlockScheduler::deliverBlock(NULL), thishandle));
  input.state = REQUESTED;
}

void BlockScheduler::deliverBlock(blkMsg *m) {
  blockMessages.push_back(m);
  CkIndex2D &src = m->src;
  int step = std::min(src.x, src.y);
  CkAssert(step == CkGetRefNum(m));
  for (StateList::iterator iter = pendingBlocks.begin(); iter != pendingBlocks.end();
       ++iter) {
    bool isU = src.x < src.y;
    BlockState::InputState &input = isU ? iter->Ustate : iter->Lstate;
    if (input.state == REQUESTED && input.m->src == src && step == iter->updatesCompleted) {
      DEBUG_SCHED("delivered %s with from (%d, %d) for (%d, %d) step = %d", isU ? "U" : "L", src.x, src.y, iter->ix, iter->iy, iter->updatesCompleted);
      input.state = ARRIVED;
      input.data = m->data;
      break;
    }
  }
  progress();
}
