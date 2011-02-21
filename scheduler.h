#ifndef LU_SCHEDULER_H
#define LU_SCHEDULER_H

#include "messages.h"
#include "luConfig.h"
#include "lu.decl.h"
#include <list>
#include <map>
#include <utility>
#include <algorithm>

enum INPUT_STATE {
  PENDING_SPACE, ALLOCATED, REQUESTED, ARRIVED
};

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
    INPUT_STATE state;
    InputState() : m(NULL), data(NULL), state(PENDING_SPACE) {}

    void reset() {
      delete m;
      m = NULL;
      data = NULL;
      state = PENDING_SPACE;
    }

    bool available() {
      return state == ARRIVED;
    }
  } Lstate, Ustate;

  BlockState(CkIndex2D index)
    : ix(index.x), iy(index.y), updatesCompleted(0), pivotsDone(false) {}
};

typedef std::list<BlockState> StateList;

class BlockScheduler : public CBase_BlockScheduler {
public:
  BlockScheduler(CProxy_LUBlk luArr_, LUConfig config)
    : luArr(luArr_), blocksHeld(0), inProgress(false) {
    blockLimit = config.memThreshold * 1024 * 1024 / (config.blockSize * config.blockSize * sizeof(double));
  }

  void registerBlock(CkIndex2D index);
  void printBlockLimit();
  void pivotsDone(CkIndex2D index);
  void dataReady(CkIndex2D index, BlockReadyMsg *m);
  void updateDone(CkIndex2D index);
  void deliverBlock(blkMsg *m);

private:
  StateList localBlocks;
  StateList pendingBlocks;
  StateList readyBlocks;
  CProxy_LUBlk luArr;
  int blockLimit, blocksHeld;
  bool inProgress;

  struct wantedBlock {
    int refs;
    bool requested;
    blkMsg* m;
    wantedBlock() : refs(0), requested(false), m(NULL) {}
  };

  std::map<std::pair<int, int>, wantedBlock> wantedBlocks;

  StateList::iterator findBlockState(CkIndex2D index);
  void progress();
  void getBlock(BlockState::InputState &input);
  bool advanceInput(BlockState::InputState &input);
  void wantBlock(BlockState::InputState &input, int x, int y);
  void incrementRefs(CkIndex2D index);

  friend class WillUse;
  friend class TryDeliver;
};

#endif
