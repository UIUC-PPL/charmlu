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
  // How many of the dependencies for the next update are not yet planned? [0..3]?
  int pendingDependencies;
  // Count of trailing updates completed and planned
  int updatesCompleted, updatesPlanned;
  // Have pivots been completed for this step?
  bool pivotsDone;
  // State of the L and U input blocks for the next update

  BlockState(CkIndex2D index)
    : ix(index.x), iy(index.y), pendingDependencies(0),
      updatesCompleted(0), updatesPlanned(0),
      pivotsDone(false)
    {}

  bool operator==(const BlockState &rhs) {
    return ix == rhs.ix && iy == rhs.iy;
  }
};

typedef std::list<BlockState> StateList;

struct Update {
  BlockState *target;
  bool triggered;
  int t;
  double *L, *U;

  void tryDeliver(int srcx, int srcy, double *data) {
    if (srcx == target->ix && srcy == t)
      L = data;
    else if (srcy == target->iy && srcx == t)
      U = data;
  }

  bool ready() { return L && U && !triggered; }

  Update(BlockState *target_, int step)
    : target(target_), triggered(false), t(step), L(NULL), U(NULL) { }

  bool operator==(const Update &rhs) {
    return target == rhs.target && t == rhs.t && L == rhs.L && U == rhs.U;
  }
};

struct ComputeU {
  int x, y, t;
  ComputeU(int x_, int y_, int t_) :
    x(x_), y(y_), t(t_) {}
};

struct PlanStep {
  int y, t;
  PlanStep(int y_, int t_) : y(y_), t(t_) {}
};

const int sdagOverheadPerBlock = 3760;

struct SchedulerProgress {
  int progress, step, allowedCols, currentPlanned;
  SchedulerProgress(int step_, int progress_, int allowedCols_,
                    int currentPlanned_)
    : step(step_), progress(progress_), allowedCols(allowedCols_),
      currentPlanned(currentPlanned_) {}
};

void registerProgressReducer();

class BlockScheduler : public CBase_BlockScheduler {
public:
  BlockScheduler(CProxy_LUBlk luArr_, LUConfig config_, CProxy_LUMgr mgr_);
  BlockScheduler(CkMigrateMessage *m) { }

  void contributeProgress(int);
  void registerBlock(CkIndex2D index);
  void allRegistered(CkReductionMsg *m);
  void incomingComputeU(CkIndex2D index, int t);
  void printBlockLimit();
  void pivotsDone(CkIndex2D index);
  void updateDone(intptr_t update_ptr);
  void factorizationDone(CkIndex2D index);
  void deliverBlock(blkMsg *m);
  void setupMulticast(rednSetupMsg *msg);
  void newColumn(CkReductionMsg *msg);

private:
  LUMgr *mgr;
  StateList localBlocks, doneBlocks;
  std::list<ComputeU> pendingComputeU;

  std::list<PlanStep> stepsToPlan;
  std::map<int, int> columnUpdatesCommitted;

  std::list<Update> plannedUpdates;
  CProxy_LUBlk luArr;
  int blockLimit;
  bool inProgress;
  int numActive;
  int totalActive;
  int previousAllowedCols;
  bool needToContribute;
  bool ownsFirstDiagonal;
  LUConfig config;
  bool countdownMode;
  int lastStep;
  int reductionCounter;

  struct wantedBlock {
    std::list<Update *> refs;
    blkMsg* m;
    wantedBlock() : m(NULL) {}
  };

  // Invariant: All of these blocks are allocated and requested
  std::map<std::pair<int, int>, wantedBlock> wantedBlocks;
  void dropRef(int srcx, int srcy, Update *update);

  std::map<std::pair<int, int>, std::list<Update*> > localWantedBlocks;

  void planUpdate(StateList::iterator target);
  void progress();
  // Makes sure that the block referenced by input is allocated
  void getBlock(int srcx, int srcy, double *&data, Update *update);
  // Potentially pick something from the ready list to update
  void runUpdate(std::list<Update>::iterator update);
};

#endif
