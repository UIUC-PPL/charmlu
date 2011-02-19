#include "lu.decl.h"
#include "lu.h"
#include <algorithm>

void BlockScheduler::registerBlock(CkIndex2D index) {
  localBlocks.push_back(BlockState(index));
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
  CkAssert(iter != localBlocks.end());
  return iter;
}

void BlockScheduler::pivotsDone(CkIndex2D index) {
  findBlockState(index)->pivotsDone = true;
}

void BlockScheduler::dataReady(CkIndex2D index, BlockReadyMsg *m) {
  StateList::iterator iter = findBlockState(index);
  BlockState::InputState &input = m->src.x == index.x ? iter->Lstate : iter->Ustate;
  input.m = m;
}

void BlockScheduler::updateDone(CkIndex2D index) {
  StateList::iterator iter = findBlockState(index);
  iter->updatesCompleted++;
  iter->pivotsDone = false;
  iter->Lstate.reset();
  iter->Ustate.reset();
}

inline bool operator==(const CkIndex2D &l, const CkIndex2D &r)
{ return l.x == r.x && l.y == r.y; }
inline bool operator<(const CkIndex2D &l, const CkIndex2D &r)
{ return l.x < r.x || (l.x == r.x && l.y < r.y); }
