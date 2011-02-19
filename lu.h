#include "lu.decl.h"
#include <vector>
#include <list>
#include <utility>
#include <map>
#include <algorithm>

struct BlockReadyMsg : public CkMcastBaseMsg, CMessage_BlockReadyMsg {
  BlockReadyMsg(CkIndex2D idx) : src(idx)
  {
    CkSetRefNum(this, std::min(src.x, src.y));
  }
  CkIndex2D src;
};

class BlockScheduler : public CBase_BlockScheduler {
public:
  BlockScheduler(CProxy_LUBlk luArr_)
    : luArr(luArr_)
  { }

  void registerBlock(CkIndex2D index);
  void pivotsDone(CkIndex2D index);
  void dataReady(CkIndex2D index, BlockReadyMsg *m);

private:
  struct BlockState {
    // Index of block
    int ix, iy;
    // Count of trailing updates completed
    int updatesCompleted;
    // Have pivots been completed for this step?
    bool pivotsDone;
    // State of the L and U input blocks for the next update
    struct InputState {
      BlockReadyMsg *m;
      double *data;
      enum {
        PENDING_SPACE, ALLOCATED, PENDING_REQUEST, REQUESTED, ARRIVED
      } state;

      InputState() : m(NULL), data(NULL), state(PENDING_SPACE) {}
    } Lstate, Ustate;

    BlockState(CkIndex2D index)
      : ix(index.x), iy(index.y), updatesCompleted(0), pivotsDone(false) {}
  };

  std::list<BlockState> localBlocks;

  CProxy_LUBlk luArr;
};
