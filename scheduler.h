#ifndef LU_SCHEDULER_H
#define LU_SCHEDULER_H

#include "messages.h"
#include "luConfig.h"
#include "lu.decl.h"
#include <list>
#include <map>
#include <utility>
#include <algorithm>

// MOVE TO scheduler.C
enum INPUT_STATE {
  PLANNED, ARRIVED, LOCAL_PLANNED
};

struct BlockState {
  // Index of block
  int ix, iy;
  // How many of the dependencies for the next update are not yet planned? 0, 1, or 2?
  int pendingDependencies;
  // Count of trailing updates completed and planned
  int updatesCompleted, updatesPlanned, updatesEligible;
  // Have pivots been completed for this step?
  bool pivotsDone;
  // State of the L and U input blocks for the next update

  BlockState(CkIndex2D index)
    : ix(index.x), iy(index.y), pendingDependencies(0),
      updatesCompleted(0), updatesPlanned(0), updatesEligible(0),
      pivotsDone(false)
    {}

  bool operator==(const BlockState &rhs) {
    return ix == rhs.ix && iy == rhs.iy;
  }
};

typedef std::list<BlockState> StateList;

//PUPbytes(std::list<Update>::iterator)

struct Update {
  BlockState *target;
  int t;
  double *L, *U;

  void tryDeliver(int srcx, int srcy, double *data) {
    if (srcx == target->ix && srcy == t) {
      CkPrintf("(%d S): Delivering L (%d,%d) to (%d,%d)\n",
	       CkMyPe(), srcx, srcy, target->ix, target->iy);
      L = data;
    }
    else if (srcy == target->iy && srcx == t) {
      CkPrintf("(%d S): Delivering U (%d,%d) to (%d,%d)\n",
	       CkMyPe(), srcx, srcy, target->ix, target->iy);
      U = data;
    } else {
      CkPrintf("(%d S): No delivery of (%d,%d) to (%d,%d)\n",
	       CkMyPe(), srcx, srcy, target->ix, target->iy);
    }
  }

  bool ready() { return L && U; }

  Update(BlockState *target_, int step)
    : target(target_), t(step), L(NULL), U(NULL) { }

  bool operator==(const Update &rhs) {
    return target == rhs.target && t == rhs.t && L == rhs.L && U == rhs.U;
  }
};

struct Panel {
  std::list<StateList::iterator> dependents;
  int updatesLeftToPlan;

  void addDependent(StateList::iterator block) { dependents.push_back(block); }
  Panel() : updatesLeftToPlan(0) { }
};

class BlockScheduler : public CBase_BlockScheduler {
public:
  BlockScheduler(CProxy_LUBlk luArr_, LUConfig config)
    : luArr(luArr_), inProgress(false) {
    blockLimit = config.memThreshold * 1024 * 1024 / (config.blockSize * config.blockSize * sizeof(double));
  }

  void registerBlock(CkIndex2D index);
  void allRegistered(CkReductionMsg *m) {
    delete m;
    progress();
  }
  void printBlockLimit();
  void pivotsDone(CkIndex2D index);
//  void dataReady(CkIndex2D index, BlockReadyMsg *m);
  void updateDone(intptr_t update_ptr);
  void factorizationDone(CkIndex2D index);
  void deliverBlock(blkMsg *m);

private:
  StateList localBlocks, doneBlocks;
  std::map<int, Panel> Lpanels;
  std::map<std::pair<int, int>, Panel> Ublocks;
  std::list<Update> plannedUpdates;
  CProxy_LUBlk luArr;
  int blockLimit;
  bool inProgress;

  struct wantedBlock {
    std::list<Update *> refs;
    blkMsg* m;
    wantedBlock() : m(NULL) {}
  };

  // Invariant: All of these blocks are allocated and requested
  std::map<std::pair<int, int>, wantedBlock> wantedBlocks;
  void dropRef(int srcx, int srcy, Update *update);

  std::map<std::pair<int, int>, std::list<Update*> > localWantedBlocks;

  void generatePlan();
  void planUpdate(StateList::iterator target);
  void progress();
  // Makes sure that the block referenced by input is allocated
  void getBlock(int srcx, int srcy, double *&data, Update *update);
  // Potentially pick something from the ready list to update
  void runUpdate(std::list<Update>::iterator update);

  template <typename K>
  void updatePanel(std::map<K, Panel> &panels, K index);
  template <typename K>
  void addDependence(std::map<K, Panel> &panels, K index, StateList::iterator block);

  // Keep localBlocks partitioned by whether blocks are eligible or not
  // Eligible at the head of localBlocks, all else at the tail
  void repositionBlock(StateList::iterator block);

#if 0
  StateList::iterator findBlockState(CkIndex2D index);
  bool advanceInput(BlockState::InputState &input, int srcX, int srcY);
  void wantBlock(BlockState::InputState &input, int x, int y);
  void incrementRefs(CkIndex2D index);
  void releaseBlock(BlockState::InputState &input);
#endif

  friend class WillUse;
  friend class TryDeliver;
};

#endif
