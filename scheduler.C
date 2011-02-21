#include "scheduler.h"
#include "lu.decl.h"
#include "messages.h"
#include "lu.h"
#include <algorithm>
#include <utility>
using std::pair;
using std::make_pair;

inline bool operator==(const CkIndex2D &l, const CkIndex2D &r)
{ return l.x == r.x && l.y == r.y; }
inline bool operator<(const CkIndex2D &l, const CkIndex2D &r)
{ return l.x < r.x || (l.x == r.x && l.y < r.y); }
pair<int, int> make_pair(CkIndex2D index) {
  return make_pair(index.x, index.y);
}

void BlockScheduler::printBlockLimit() {
  CkPrintf("%d: block limit = %d\n", CkMyPe(), blockLimit);
}

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
  DEBUG_SCHED("dataReady message src = (%d, %d)", input.m->src.x, input.m->src.y);
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

struct earliestRelevantSorter {
  bool operator()(const BlockState& state1, const BlockState& state2) {
    return state1.ix != state2.ix ? state1.ix < state2.ix : state1.iy < state2.iy;
  }
};

bool BlockScheduler::advanceInput(BlockState::InputState &input) {
    bool stateModified = false;

    if (input.state == PENDING_SPACE) {
	input.state = ALLOCATED;
	stateModified = true;
    }

    if (input.state == ALLOCATED) {
	if (input.m) {
            getBlock(input);
            stateModified = true;
	} // else we wait for ready msg
    }

    return stateModified;
}

void BlockScheduler::wantBlock(BlockState::InputState &input, int x, int y) {
  wantedBlocks.insert(make_pair(make_pair(x, y), wantedBlock()));
}

void BlockScheduler::progress() {
  // Prevent reentrance
  if (inProgress)
    return;

  inProgress = true;
  bool stateModified;

  do {
    stateModified = false;
    localBlocks.sort(updatesCompletedSorter());

    while (wantedBlocks.size() < blockLimit && localBlocks.size() > 0) {
      // Move some blocks from localBlocks to pendingBlocks
      BlockState &block = *localBlocks.begin();
      DEBUG_SCHED("putting into pendingBlocks: (%d, %d)", block.ix, block.iy);
      wantBlock(block.Lstate, block.ix, block.updatesCompleted);
      wantBlock(block.Ustate, block.updatesCompleted, block.iy);
      pendingBlocks.splice(pendingBlocks.end(), localBlocks, localBlocks.begin());
      stateModified = true;
    }

    // Try to advance state for each
    for (StateList::iterator iter = pendingBlocks.begin(); iter != pendingBlocks.end();
         ++iter) {
      BlockState &block = *iter;
      DEBUG_SCHED("examining pending block (%d, %d), Lstate = %d, Ustate = %d", block.ix, block.iy,
                  block.Lstate.state, block.Ustate.state);
      if (block.pivotsDone) {
	stateModified |= advanceInput(block.Lstate);
	stateModified |= advanceInput(block.Ustate);

        if (block.Ustate.available() && block.Lstate.available()) {
          DEBUG_SCHED("both block indicate arrived for (%d, %d)", block.ix, block.iy);
          readyBlocks.push_back(block);
          iter = pendingBlocks.erase(iter);
          stateModified = true;
        }
      }
    }

    StateList::iterator block = readyBlocks.begin();
    if (block != readyBlocks.end()) {
      CkIndex2D Lsrc = block->Lstate.m->src, Usrc = block->Ustate.m->src;
      CkAssert(block->updatesCompleted == Lsrc.y && block->updatesCompleted == Usrc.x);

      int Lref = --wantedBlocks[make_pair(block->ix, block->updatesCompleted)].refs;
      int Uref = --wantedBlocks[make_pair(block->updatesCompleted, block->iy)].refs;

      CkAssert(block->Lstate.data && block->Ustate.data);

      luArr(block->ix, block->iy).ckLocal()->
        processTrailingUpdate(block->updatesCompleted, (intptr_t)block->Lstate.data,
                              (intptr_t)block->Ustate.data);

      if (Lref == 0) {
        std::map<pair<int, int>, wantedBlock>::iterator iter =
          wantedBlocks.find(make_pair(block->ix, block->updatesCompleted-1));
        delete iter->second.m;
        wantedBlocks.erase(iter);
      }
      if (Uref == 0) {
        std::map<pair<int, int>, wantedBlock>::iterator iter =
          wantedBlocks.find(make_pair(block->updatesCompleted-1, block->iy));
        delete iter->second.m;
        wantedBlocks.erase(iter);
      }

      if (block->updatesCompleted < std::min(block->ix, block->iy))
        localBlocks.splice(localBlocks.end(), readyBlocks, block);
      else
        readyBlocks.erase(block);

      stateModified = true;
    }
  } while (stateModified);

  inProgress = false;
}

struct WillUse {
  BlockScheduler::wantedBlock &block;
  pair<int, int> index;
  WillUse(BlockScheduler::wantedBlock &block_, pair<int, int> index_) : block(block_), index(index_) {}
  void operator()(BlockState &rblock) {
    if (rblock.ix == index.first && index.second > rblock.iy ||
        rblock.iy == index.second && index.first > rblock.ix) {
      block.refs++;
    }
  }
};

void BlockScheduler::getBlock(BlockState::InputState &input) {
  CkAssert(input.state == ALLOCATED);

  pair<int, int> src = make_pair(input.m->src);

  input.state = REQUESTED;

  if (wantedBlocks[src].m) {
    DEBUG_SCHED("already ARRIVED from (%d, %d)", input.m->src.x, input.m->src.y);
    input.data = wantedBlocks[src].m->data;
    input.state = ARRIVED;
    return;
  }
  if (wantedBlocks[src].requested) {
    DEBUG_SCHED("already REQUESTED from (%d, %d)", input.m->src.x, input.m->src.y);
    return;
  }

  std::for_each(pendingBlocks.begin(), pendingBlocks.end(), WillUse(wantedBlocks[src], src));
  std::for_each(localBlocks.begin(), localBlocks.end(), WillUse(wantedBlocks[src], src));
  std::for_each(readyBlocks.begin(), readyBlocks.end(), WillUse(wantedBlocks[src], src));

  wantedBlocks[src].requested = true;
  DEBUG_SCHED("requesting getBlock from (%d, %d)", input.m->src.x, input.m->src.y);
  luArr(src.first, src.second).getBlock(CkCallback(CkIndex_BlockScheduler::deliverBlock(NULL), thishandle));
}

struct TryDeliver {
  BlockScheduler::wantedBlock &block;
  TryDeliver(BlockScheduler::wantedBlock &block_) : block(block_) {}
  void operator()(BlockState &rblock) {
    const CkIndex2D src = block.m->src;
    int step = std::min(src.x, src.y);
    bool isU = src.x < src.y;
    BlockState::InputState &input = isU ? rblock.Ustate : rblock.Lstate;
    if (input.m && input.m->src == src && step == rblock.updatesCompleted) {
      DEBUG_SCHED("delivered %s with from (%d, %d) for (%d, %d) step = %d", isU ? "U" : "L",
                  src.x, src.y, rblock.ix, rblock.iy, rblock.updatesCompleted);
      input.state = ARRIVED;
      input.data = block.m->data;
    }
  }
};

void BlockScheduler::deliverBlock(blkMsg *m) {
  pair<int, int> src = make_pair(m->src);
  wantedBlocks[src].m = m;
  std::for_each(pendingBlocks.begin(), pendingBlocks.end(), TryDeliver(wantedBlocks[src]));
  std::for_each(localBlocks.begin(), localBlocks.end(), TryDeliver(wantedBlocks[src]));
  progress();
}
