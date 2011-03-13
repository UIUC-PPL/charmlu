#include "scheduler.h"
#include "lu.decl.h"
#include "messages.h"
#include "lu.h"
#include <algorithm>
using std::min;
#include <utility>
using std::pair;
using std::make_pair;

inline bool operator==(const CkIndex2D &l, const CkIndex2D &r)
{ return l.x == r.x && l.y == r.y; }
pair<int, int> make_pair(CkIndex2D index) {
  return make_pair(index.x, index.y);
}

BlockScheduler::BlockScheduler(CProxy_LUBlk luArr_, LUConfig config, CProxy_LUMgr mgr_)
  : luArr(luArr_), mgr(mgr_.ckLocalBranch()), inProgress(false), numActive(0),
    pendingTriggered(0), sendDelay(0) {
  blockLimit = config.memThreshold * 1024 * 1024 /
    (config.blockSize * (config.blockSize + 1) * sizeof(double) + sizeof(LUBlk) + sdagOverheadPerBlock);

  contribute(CkCallback(CkIndex_LUBlk::schedulerReady(NULL), luArr));
}

void BlockScheduler::incomingComputeU(CkIndex2D index, int t) {
  if (numActive > 0) {
    pendingComputeU.push_back(ComputeU(index.x, index.y, t));
  } else {
    CkEntryOptions opts;
    luArr(index).processComputeU(0, &(mgr->setPrio(RECVL, opts, index.y)));
    pendingTriggered++;
  }
}

void BlockScheduler::scheduleSend(CkIndex2D sender) {
  std::list<CkIndex2D>::iterator iter = find(scheduledSends.begin(),
                                             scheduledSends.end(), sender);

  if (iter == scheduledSends.end()) {
    if (pendingTriggered == 0) {
      CkEntryOptions opts;
      luArr[sender].sendBlocks(0, &mgr->setPrio(SEND_BLOCKS, opts));
    } else {
      scheduledSends.push_back(sender);
    }
  }
}

void BlockScheduler::printBlockLimit() {
  CkPrintf("%d: block limit = %d\n", CkMyPe(), blockLimit);
}

void BlockScheduler::registerBlock(CkIndex2D index) {
  blockLimit--;
  if (index.x != 0 && index.y != 0) {
    localBlocks.push_back(BlockState(index));
    for (int i = 1; i < min(index.x, index.y) + 1; i++) {
      panels[i].updatesLeftToPlan++;
    }
  }

  for (int i = 0; i < index.x; i++)
    luArr(i, index.y).prepareForMulticast(CkMyPe());
  for (int i = 0; i < index.y; i++)
    luArr(index.x, i).prepareForMulticast(CkMyPe());

  if (blockLimit< 2)
    CkAbort("Too little space to plan even one trailing update");
}

void BlockScheduler::allRegistered(CkReductionMsg *m) {
  delete m;
  progress();
}

void BlockScheduler::startedActivePanel() {
  numActive++;
}

void BlockScheduler::setupMulticast(rednSetupMsg *msg) {
  CkMulticastMgr* mcastMgr = CProxy_CkMulticastMgr(msg->rednMgrGID).ckLocalBranch();
  CkSectionInfo cookie;
  CkGetSectionInfo(cookie, msg);
  mcastMgr->contribute(0, NULL, CkReduction::sum_int, cookie);
  delete msg;
}

void BlockScheduler::repositionBlock(StateList::iterator block) {
  StateList::iterator pos = block->pendingDependencies == 0 ?
    localBlocks.begin() : localBlocks.end();
  localBlocks.splice(pos, localBlocks, block);
}

template <typename K>
void BlockScheduler::updatePanel(std::map<K, Panel> &panels, K index) {
  typename std::map<K, Panel>::iterator iter = panels.find(index);
  CkAssert(iter != panels.end());
  Panel &panel = iter->second;

  panel.updatesLeftToPlan--;
  if(panel.updatesLeftToPlan == 0) {
    for (std::list<StateList::iterator>::iterator i = panel.dependents.begin();
         i != panel.dependents.end(); ++i) {
      (*i)->pendingDependencies--;
      repositionBlock(*i);
    }

    panels.erase(iter);
  }
}

template <typename K>
void BlockScheduler::addDependence(std::map<K, Panel> &panels, K index,
				   StateList::iterator block) {
  typename std::map<K, Panel>::iterator panel = panels.find(index);
  if (panel != panels.end()) {
    block->pendingDependencies++;
    panel->second.addDependent(block);
  }
}

void BlockScheduler::planUpdate(StateList::iterator target) {
  CkAssert(target->pendingDependencies == 0);

  int t = target->updatesPlanned++;

  plannedUpdates.push_back(Update(&*target, t));
  Update &update = plannedUpdates.back();

  getBlock(target->ix, t, update.L, &update);
  getBlock(t, target->iy, update.U, &update);

  if (target->updatesPlanned < min(target->ix, target->iy)) {
    addDependence(panels, t+1, target);
    repositionBlock(target);
  } else {
    doneBlocks.splice(doneBlocks.end(), localBlocks, target);
  }

  updatePanel(panels, t+1);
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
    update->tryDeliver(srcx, srcy, block.data);
  }
}

void BlockScheduler::deliverBlock(blkMsg *m) {
  DEBUG_SCHED("deliverBlock src (%d, %d)", m->src.x, m->src.y);

  wantedBlock &block = wantedBlocks[make_pair(m->src)];
  block.m = m;
  block.data = m->data;

  for (std::list<Update*>::iterator update = block.refs.begin();
       update != block.refs.end(); ++update) {
    DEBUG_SCHED("tryDeliver of (%d,%d) to (%d,%d) @ %d", m->src.x, m->src.y,
		(*update)->target->ix, (*update)->target->iy, (*update)->t);
    (*update)->tryDeliver(m->src.x, m->src.y, m->data);
  }

  CkAssert(m->npes >= 1);
  CkAssert(CkMyPe() == m->pes[m->offset]);

  // This processor is no longer part of the set
  m->offset++;
  m->npes--;
  propagateBlkMsg(m, thisProxy);

  progress();
}

void propagateBlkMsg(blkMsg *m, CProxy_BlockScheduler bs) {
  DEBUG_SCHED("Delivering to processors");
  for (int i = m->offset; i < m->npes; ++i)
    DEBUG_SCHED("\t%d", m->pes[i]);

  if (m->npes >= 2) {
    blkMsg *m2 = (blkMsg *)CkCopyMsg((void **)&m);
    m2->offset += m->npes / 2;
    m2->npes -= m->npes / 2;
    m->npes -= m2->npes;

    bs[m2->pes[m2->offset]].deliverBlock(m2);
  }

  if (m->npes >= 1) {
    takeRef(m);
    bs[m->pes[m->offset]].deliverBlock(m);
  }

  traceMemoryUsage();
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
  pendingTriggered++;

  CkEntryOptions opts;
  int t = update.t;
  intptr_t update_ptr = (intptr_t)&update;
  luArr(tx, ty).processTrailingUpdate(t, update_ptr,
                                      &(mgr->setPrio(PROCESS_TRAILING_UPDATE, opts, ty, tx)));
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

  pendingTriggered--;
  CkAssert(pendingTriggered >= 0);

  plannedUpdates.erase(std::find(plannedUpdates.begin(), plannedUpdates.end(), update));

  runScheduledSends();
  progress();
}

void BlockScheduler::updateUntriggered() {
  pendingTriggered--;
  if (pendingTriggered == 0) {
    runScheduledSends();
  }
}

static const int SEND_SKIP = 10;

bool stepColumnOrder(const CkIndex2D &l, const CkIndex2D &r) {
  return min(l.x, l.y) != min(r.x, r.y) ? min(l.x,l.y) < min(r.x, r.y) : l.y < r.y;
}

void BlockScheduler::runScheduledSends() {
  if (scheduledSends.size() > 0) {
    sendDelay++;
    if (pendingTriggered != 0 && sendDelay >= SEND_SKIP) {
      CkEntryOptions opts;
      scheduledSends.sort(stepColumnOrder);
      luArr[scheduledSends.front()].sendBlocks(0, &mgr->setPrio(SEND_BLOCKS, opts));
      scheduledSends.pop_front();
      sendDelay = 0;
    } else if (pendingTriggered == 0) {
      for (std::list<CkIndex2D>::iterator iter = scheduledSends.begin();
       iter != scheduledSends.end(); ++iter) {
        CkEntryOptions opts;
        luArr[*iter].sendBlocks(0, &mgr->setPrio(SEND_BLOCKS, opts));
      }
      scheduledSends.clear();
    }
  }
}

bool BlockScheduler::shouldExecute() {
  return numActive == 0;
}

void BlockScheduler::factorizationDone(CkIndex2D index) {
  DEBUG_SCHED("factorizationDone on (%d,%d)", index.x, index.y);

  if (index.x >= index.y)
    numActive--;

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

void BlockScheduler::progress() {
  DEBUG_SCHED("Called progress, already? %s", inProgress? "true" : "false" );
  // Prevent reentrance
  if (inProgress)
    return;

  inProgress = true;

  bool stateModified;

  do {
    stateModified = false;
    bool plannedAnything = true;
    while (wantedBlocks.size() < blockLimit && plannedAnything) {
      plannedAnything = false;
      if (localBlocks.size() > 0) {
	DEBUG_SCHED("Local Blocks: ");
	for (StateList::iterator block = localBlocks.begin(); block != localBlocks.end();
	     ++block) {
	  DEBUG_SCHED("\t(%d,%d), deps: %d, updatesCompleted %d, updatesPlanned %d",
		      block->ix, block->iy, block->pendingDependencies,
		      block->updatesCompleted, block->updatesPlanned);
	}

        localBlocks.sort(eligibilityYOrder);

	CkAssert(localBlocks.front().pendingDependencies == 0);
	planUpdate(localBlocks.begin());
	plannedAnything = true;
      }
    }

    if (numActive == 0) {
      // Start processComputeU
      // TODO: refactor into a foreach?
      for (std::list<ComputeU>::iterator computeU = pendingComputeU.begin();
           computeU != pendingComputeU.end(); ++computeU) {
        CkEntryOptions opts;
        luArr(computeU->x, computeU->y).
          processComputeU(0, &(mgr->setPrio(RECVL, opts, computeU->y)));
        pendingTriggered++;
        computeU = pendingComputeU.erase(computeU);
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
    }
  } while (stateModified);

  inProgress = false;
}

