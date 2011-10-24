#ifndef LU_SCHEDULER_H
#define LU_SCHEDULER_H

#include "messages.h"
#include "luConfig.h"
#include "lu.decl.h"
#include <list>
#include <map>
#include <utility>
#include <algorithm>
#include <utility>

enum INPUT_STATE {
  PLANNED, ARRIVED, LOCAL_PLANNED
};

class Panel;

struct BlockState {
  // Index of block
  int ix, iy;
  // How many of the dependencies for the next update are not yet planned? [0..3]?
  std::list<Panel*> pendingDependencies;
  // Count of trailing updates completed and planned
  int updatesCompleted, updatesPlanned, updatesEligible;
  // Have pivots been completed for this step?
  bool pivotsDone;

  BlockState(CkIndex2D index)
    : ix(index.x), iy(index.y), pendingDependencies(0),
      updatesCompleted(0), updatesPlanned(0), updatesEligible(0),
      pivotsDone(false)
  { }

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

  bool isComputeU() { return target->ix == t; }
  bool ready() { return L && (U || isComputeU()) && !triggered; }

  Update(BlockState *target_, int step)
    : target(target_), triggered(false), t(step), L(NULL), U(NULL) { }

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

struct ComputeU {
  int x, y, t;
  ComputeU(int x_, int y_, int t_) :
    x(x_), y(y_), t(t_) {}
};

const int sdagOverheadPerBlock = 3760;

class BlockScheduler : public CBase_BlockScheduler {
public:
  BlockScheduler(CProxy_LUBlk luArr_, LUConfig config, CProxy_LUMgr mgr_);
  BlockScheduler(CkMigrateMessage *m) { }

  void registerBlock(CkIndex2D index);
  void allRegistered(CkReductionMsg *m);
  void printBlockLimit();
  void pivotsDone(CkIndex2D index);
  void updateDone(intptr_t update_ptr);
  void factorizationDone(CkIndex2D index);
  void deliverBlock(blkMsg *m);
  void startedActivePanel();
  bool shouldExecute();
  void scheduleSend(blkMsg *m, bool onActive);
  void scheduleSend(CkIndex2D index, bool onActive);
  void updateUntriggered();
  void pumpMessages();
  void releaseActiveColumn(const int y, const int t);
  void outputStats();

private:
  LUMgr *mgr;
  StateList localBlocks, doneBlocks;

  std::map<int, Panel> panels, Upanels;
  std::list<Update> plannedUpdates;
  CProxy_LUBlk luArr;
  int blockLimit;
  bool inProgress, inPumpMessages;
  int numActive;
  int pendingTriggered;
  bool reverseSends;

  int sendDelay;
  std::list<blkMsg *> scheduledSends, sendsInFlight;

  struct wantedBlock {
    std::list<Update *> refs;
    blkMsg* m;
    double *data;
    bool isSending;
    wantedBlock() : m(NULL), data(NULL), isSending(false) {}
  };

  void propagateBlkMsg(blkMsg *msg);
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

  template <typename K>
  void updatePanel(std::map<K, Panel> &panels, K index);
  template <typename K>
  void addDependence(std::map<K, Panel> &panels, K index, StateList::iterator block);

  // Memory usage instrumentation
  // Memory occupied before factorization starts (kilobytes)
  int baseMemory;
  // Maximum memory usage seen (kilobytes)
  int maxMemory, maxMemoryIncreases, maxMemoryStep;
};

#endif
