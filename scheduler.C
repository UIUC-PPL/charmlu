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

BlockScheduler::BlockScheduler(CProxy_LUBlk luArr_, LUConfig config_, CProxy_LUMgr mgr_)
  : luArr(luArr_), mgr(mgr_.ckLocalBranch()), inProgress(false), config(config_),
    numActive(0), ownsFirstDiagonal(false), previousAllowedCols(1),
    needToContribute(false), countdownMode(false), lastStep(-1) {
  blockLimit = config.memThreshold * 1024 * 1024 /
    (config.blockSize * (config.blockSize + 1) * sizeof(double) + sizeof(LUBlk) + sdagOverheadPerBlock);

  contribute(CkCallback(CkIndex_LUBlk::schedulerReady(NULL), luArr));
}

void BlockScheduler::incomingComputeU(CkIndex2D index, int t) {
  // std::map<int, int>::iterator apanel = activePanels.find(t + 1);
  // if (apanel != activePanels.end() && apanel->second > 0 &&
  //     index.y != t + 1) {
  //   pendingComputeU.push_back(ComputeU(index.x, index.y, t));
  // } else {
  //   CkEntryOptions opts;
  //   luArr(index).processComputeU(0, &(mgr->setPrio(RECVL, opts, index.y)));
  // }
}

void BlockScheduler::printBlockLimit() {
  CkPrintf("%d: block limit = %d\n", CkMyPe(), blockLimit);
}

void BlockScheduler::registerBlock(CkIndex2D index) {
  blockLimit--;
  if (index.x != 0 && index.y != 0) {
    localBlocks.push_back(BlockState(index));
  }

  if (index.x == index.y && index.x == 0)
    ownsFirstDiagonal = true;

  if (index.y != 0)
    columnUpdatesCommitted[index.y] = -1;

  for (int i = 0; i < index.x; i++)
    luArr(i, index.y).prepareForMulticast(CkMyPe());
  for (int i = 0; i < index.y; i++)
    luArr(index.x, i).prepareForMulticast(CkMyPe());

  if (blockLimit< 2)
    CkAbort("Too little space to plan even one trailing update");
}

void BlockScheduler::allRegistered(CkReductionMsg *m) {
  delete m;
  contributeProgress(0);
}

void BlockScheduler::setupMulticast(rednSetupMsg *msg) {
  CkMulticastMgr* mcastMgr = CProxy_CkMulticastMgr(msg->rednMgrGID).ckLocalBranch();
  CkSectionInfo cookie;
  CkGetSectionInfo(cookie, msg);
  mcastMgr->contribute(0, NULL, CkReduction::sum_int, cookie);
  delete msg;
}

void BlockScheduler::planUpdate(StateList::iterator target) {
  int t = target->updatesPlanned++;

  plannedUpdates.push_back(Update(&*target, t));
  Update &update = plannedUpdates.back();

  getBlock(target->ix, t, update.L, &update);
  getBlock(t, target->iy, update.U, &update);

  CkAssert(target->updatesPlanned <= min(target->ix, target->iy));

  if (target->updatesPlanned == min(target->ix, target->iy))
    doneBlocks.splice(doneBlocks.end(), localBlocks, target);
}

void BlockScheduler::getBlock(int srcx, int srcy, double *&data,
			      Update *update) {
  LUBlk *local = luArr(srcx, srcy).ckLocal();
  if (local) {
    if (local->factored) {
      DEBUG_SCHED("Found block (%d, %d) ready locally for update to (%d,%d)",
		  srcx, srcy, update->target->ix, update->target->iy);
      data = local->getBlock();
    } else {
      DEBUG_SCHED("Found block (%d, %d) pending locally for update to (%d,%d)",
		  srcx, srcy, update->target->ix, update->target->iy);
      localWantedBlocks[make_pair(srcx, srcy)].push_back(update);
    }
    return;
  }

  pair<int, int> src = make_pair(srcx, srcy);
  wantedBlock &block = wantedBlocks[src];

  block.refs.push_back(update);

  if (block.refs.size() == 1) {
    // First reference to this block, so ask for it
    DEBUG_SCHED("requesting getBlock from (%d, %d)", srcx, srcy);
    CkEntryOptions opts;
    luArr(src.first, src.second).
      getBlock(CkMyPe(), &(mgr->setPrio(GET_BLOCK, opts)));
  }

  if (block.m) {
    DEBUG_SCHED("already ARRIVED from (%d, %d)", srcx, srcy);
    update->tryDeliver(srcx, srcy, block.m->data);
  }
}

void BlockScheduler::deliverBlock(blkMsg *m) {
  DEBUG_SCHED("deliverBlock src (%d, %d)", m->src.x, m->src.y);
  wantedBlock &block = wantedBlocks[make_pair(m->src)];
  block.m = m;

  for (std::list<Update*>::iterator update = block.refs.begin();
       update != block.refs.end(); ++update) {
    DEBUG_SCHED("tryDeliver of (%d,%d) to (%d,%d) @ %d", m->src.x, m->src.y,
		(*update)->target->ix, (*update)->target->iy, (*update)->t);
    (*update)->tryDeliver(m->src.x, m->src.y, m->data);
  }

  progress();
}

void BlockScheduler::dropRef(int srcx, int srcy, Update *update) {
  std::map<std::pair<int, int>, wantedBlock>::iterator input =
    wantedBlocks.find(make_pair(srcx, srcy));

  if (input == wantedBlocks.end())
    return;

  input->second.refs.remove(update);
  if (input->second.refs.size() == 0) {
    delete input->second.m;
    wantedBlocks.erase(input);
  }
}

void BlockScheduler::runUpdate(std::list<Update>::iterator iter) {
  Update &update = *iter;
  CkAssert(update.ready());
  int tx = update.target->ix, ty = update.target->iy;
  update.triggered = true;

  CkEntryOptions opts;
  int t = update.t;
  intptr_t update_ptr = (intptr_t)&update;
  luArr(tx, ty).processTrailingUpdate(t, update_ptr,
                                      &(mgr->setPrio(PROCESS_TRAILING_UPDATE, opts, ty)));
}

void BlockScheduler::updateDone(intptr_t update_ptr) {
  Update &update = *(Update *)update_ptr;
  int tx = update.target->ix, ty = update.target->iy;
  DEBUG_SCHED("updateDone on (%d,%d)", tx, ty);

  dropRef(tx, update.t, &update);
  dropRef(update.t, ty, &update);

  update.target->updatesCompleted++;
  if (update.target->updatesCompleted == min(update.target->ix, update.target->iy)) {
    // Last update on this block
    //doneBlocks.erase(std::find(doneBlocks.begin(), doneBlocks.end(), *update.target));
  }

  plannedUpdates.erase(std::find(plannedUpdates.begin(), plannedUpdates.end(), update));

  if (countdownMode)
    contributeProgress(0);

  progress();
}

void BlockScheduler::factorizationDone(CkIndex2D index) {
  DEBUG_SCHED("factorizationDone on (%d,%d)", index.x, index.y);

  if (index.x == index.y && needToContribute) {
    contributeProgress(0);
    needToContribute = false;
  }

  std::map<std::pair<int, int>, std::list<Update*> >::iterator wanters =
    localWantedBlocks.find(make_pair(index.x, index.y));
  if (wanters != localWantedBlocks.end()) {
    std::list<Update*> &wantList = wanters->second;
    CkAssert(luArr(index).ckLocal());

    for (std::list<Update*>::iterator wanter = wantList.begin();
	 wanter != wantList.end(); ++wanter) {
      DEBUG_SCHED("tryDeliver of (%d,%d) to (%d,%d) @ %d", index.x, index.y,
		  (*wanter)->target->ix, (*wanter)->target->iy, (*wanter)->t);
      (*wanter)->tryDeliver(index.x, index.y, luArr(index).ckLocal()->getBlock());
    }

    localWantedBlocks.erase(wanters);
  }

  progress();
}

bool eligibilityYOrder(const BlockState& block1, const BlockState& block2) {
  if (block1.pendingDependencies != block2.pendingDependencies) {
    return block1.pendingDependencies < block2.pendingDependencies;
  } else {
    return block1.iy < block2.iy;
  }
}

const int PIPE_SIZE = 2;

void BlockScheduler::newColumn(CkReductionMsg *msg) {
  SchedulerProgress sp = *(SchedulerProgress*)msg->getData();

  DEBUG_SCHED("newColumn current active panel = %d, reduction num = %d, progress = %d",
              sp.step, msg->getRedNo(), sp.progress);

  if (sp.step == config.numBlocks - 1) {
    CkPrintf("Shutting DOWN, redNo = %d\n", msg->getRedNo());
    return;
  }

  std::map<int, int>::iterator iter;

  if (lastStep != sp.step || sp.currentPlanned == 0) {
    for (iter = columnUpdatesCommitted.begin();
         iter != columnUpdatesCommitted.end(); ++iter) {
      DEBUG_SCHED("Considering column (y = %d, t = %d)", iter->first, iter->second);
      if (iter->second < sp.step) {
        DEBUG_SCHED("Determined new column to plan: (y = %d, t = %d), reduction num = %d",
                    iter->first, iter->second + 1, msg->getRedNo());
        iter->second++;
        stepsToPlan.push_back(PlanStep(iter->first, iter->second));
        break;
      }
    }
  }

  delete msg;

  int allowedCols = 0;

  for (iter = columnUpdatesCommitted.begin();
       iter != columnUpdatesCommitted.end(); ++iter) {
    if (iter->second < sp.step) {
      allowedCols++;
    }
  }

  previousAllowedCols = allowedCols;

  DEBUG_SCHED("lastStep = %d, sp.step = %d", lastStep, sp.step);
#if 0
  if (lastStep != sp.step) {
    /*if (lastStep != sp.step) {
      reductionCounter = 2;
    } else {
      reductionCounter--;
      }*/
    if (sp.allowedCols != 0) {
      contributeProgress(0);
    } else { 
      if (luArr(sp.step+1, sp.step+1).ckLocal() &&
          !luArr(sp.step+1, sp.step+1).ckLocal()->factored) {
        needToContribute = true;
      } else {
        contributeProgress(0);
      }
    }
  } else {
  }
#endif

  bool foundTriggered = false;
  for (std::list<Update>::iterator update = plannedUpdates.begin();
       update != plannedUpdates.end(); ++update) {
    if (update->triggered)
      foundTriggered = true;
  }

  if (foundTriggered) {
    countdownMode = true;
    //needToContribute = true;
  } else {
    if (sp.allowedCols != 0) {
      contributeProgress(0);
    } else { 
      if (luArr(sp.step+1, sp.step+1).ckLocal() &&
          !luArr(sp.step+1, sp.step+1).ckLocal()->factored) {
        needToContribute = true;
      } else {
        contributeProgress(0);
      }
    }
  }

  lastStep = sp.step;

  progress();
}

static CkReductionMsg* progressReducerFn(int count, CkReductionMsg **msgs) {
  SchedulerProgress sp = *(SchedulerProgress*)msgs[0]->getData();
  for (int i = 1; i < count; i++) {
    SchedulerProgress *cur = (SchedulerProgress*)msgs[i]->getData();
    if (cur->step > sp.step) {
      sp.step = cur->step;
      sp.progress = cur->progress;
    }
    sp.allowedCols += cur->allowedCols;
    sp.currentPlanned += cur->currentPlanned;
  }
  DEBUG_SCHED("progressReducerFn returning: step = %d, progress = %d, allowedCols = %d",
              sp.step, sp.progress, sp.allowedCols);
  return CkReductionMsg::buildNew(sizeof(SchedulerProgress), &sp);
}

static CkReduction::reducerType progressReducer;

void registerProgressReducer() {
  progressReducer = CkReduction::addReducer(progressReducerFn);
}

void BlockScheduler::contributeProgress(int) {
  int progress = 0, step = 0;
  for (StateList::iterator block = doneBlocks.begin();
       block != doneBlocks.end(); ++block) {
    if (block->ix == block->iy && block->updatesCompleted == block->ix &&
        block->updatesCompleted > step) {
      step = block->updatesCompleted;
      progress = 0;
      //progress = luArr(block->ix, block->iy).ckLocal()->getActivePanelProgress();
    }
  }

  DEBUG_SCHED("contributing scheduler progress: step = %d, progress = %d", step, progress);

  SchedulerProgress schedProgress(step, progress, previousAllowedCols, plannedUpdates.size());
  contribute(sizeof(schedProgress), &schedProgress, progressReducer,
             CkCallback(CkIndex_BlockScheduler::newColumn(NULL), thisProxy));
}

struct InChosenColumn {
  PlanStep step;
  InChosenColumn(PlanStep step_) : step(step_) {}
  bool operator()(BlockState &block) {
    return block.iy == step.y && block.updatesPlanned == step.t;
  }
};

void BlockScheduler::progress() {
  DEBUG_SCHED("Called progress, already? %s", inProgress? "true" : "false" );
  // Prevent reentrance
  if (inProgress)
    return;

  inProgress = true;

  // if (countdownContribute == 0 && countdownMode) {
  //   contributeProgress(0);
  // }

  bool stateModified;

  do {
    stateModified = false;
    bool plannedAnything = true;
    while (wantedBlocks.size() < 1 && stepsToPlan.size() > 0 && plannedAnything) {
      plannedAnything = false;
      if (localBlocks.size() > 0) {
	DEBUG_SCHED("Local Blocks: ");
	for (StateList::iterator block = localBlocks.begin(); block != localBlocks.end();
	     ++block) {
	  DEBUG_SCHED("\t(%d,%d), deps: %d, updatesCompleted %d, updatesPlanned %d",
		      block->ix, block->iy, block->pendingDependencies,
		      block->updatesCompleted, block->updatesPlanned);
	}

        PlanStep nextStep = stepsToPlan.front();

        StateList::iterator block =
          find_if(localBlocks.begin(), localBlocks.end(), InChosenColumn(nextStep));

        if (block != localBlocks.end()) {
          DEBUG_SCHED("Planning update (%d, %d) for %d, nextStep: y = %d, t = %d",
                      block->ix, block->iy, block->updatesPlanned, nextStep.y, nextStep.t);
          planUpdate(block);
        } else {
          stepsToPlan.pop_front();
          //contributeProgress(0);
        }
        plannedAnything = true;
      }
    }

    // Start trailing updates
    for (std::list<Update>::iterator update = plannedUpdates.begin();
         update != plannedUpdates.end(); ++update) {
      if (update->ready()) {
        runUpdate(update);
        stateModified = true;
        break;
      }
    }
  } while (stateModified);

  inProgress = false;
}

