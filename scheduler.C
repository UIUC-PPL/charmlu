#include "scheduler.h"
#include "lu.decl.h"
#include "messages.h"
#include "lu.h"
#include <algorithm>
using std::min;
using std::find;
#include <utility>
using std::pair;
using std::make_pair;
#include "register.h"
using std::list;
using std::map;

inline bool operator==(const CkIndex2D &l, const CkIndex2D &r)
{ return l.x == r.x && l.y == r.y; }
pair<int, int> make_pair(CkIndex2D index) {
  return make_pair(index.x, index.y);
}

BlockScheduler::BlockScheduler(CProxy_LUBlk luArr_, LUConfig config, CProxy_LUMgr mgr_)
  : luArr(luArr_), mgr(mgr_.ckLocalBranch()), inProgress(false), numActive(0),
    pendingTriggered(0), sendDelay(0), reverseSends(CkMyPe() % 2 == 0) {
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

void BlockScheduler::scheduleSend(blkMsg *msg) {
  list<blkMsg *>::iterator iter = find(scheduledSends.begin(),
                                       scheduledSends.end(), msg);
  DEBUG_SCHED("scheduling a new send (%d, %d)", msg->src.x, msg->src.y);
  if (iter == scheduledSends.end()) {
    scheduledSends.push_back(msg);
  }

  pumpMessages();
}

void pumpOnIdle(void *s, double) {
  BlockScheduler *scheduler = (BlockScheduler *)s;
  //CkPrintf("Firing pumpMessages on Idle\n");
  //CcdCallOnCondition(CcdPERIODIC_100ms, pumpOnIdle, s);
  scheduler->pumpMessages();
}

static const int SEND_LIMIT = 1;

void BlockScheduler::pumpMessages() {
  for (list<blkMsg*>::iterator iter = sendsInFlight.begin();
       iter != sendsInFlight.end(); ++iter) {
    //DEBUG_SCHED("pumpMessages, iter through sendsInFlight %p", *iter);
    blkMsg *msg = *iter;
    int ref = REFFIELD(UsrToEnv(msg));
    CkAssert(ref > 0 && ref <= 2);
    if (ref == 1) {
      if (!msg->firstHalfSent) {
        DEBUG_SCHED("%p calling propagate on second half", *iter);
        msg->firstHalfSent = true;
        propagateBlkMsg(msg);
      } else {
        map<pair<int, int>, wantedBlock>::iterator blockIter =
          wantedBlocks.find(make_pair(msg->src));
        if (blockIter != wantedBlocks.end()) {
          wantedBlock &block = blockIter->second;
          CkAssert(block.isSending);
          block.isSending = false;
          if (block.refs.size() == 0)  {
            delete msg;
            wantedBlocks.erase(blockIter);
          }
        } else {
          // If local, needs to be set back to false, so setupMsg gets called
          msg->firstHalfSent = false;
        }
        iter = sendsInFlight.erase(iter);
      }
    }
  }

  for (list<blkMsg*>::iterator iter = scheduledSends.begin();
       iter != scheduledSends.end() && sendsInFlight.size() < SEND_LIMIT;
       ++iter) {
    //DEBUG_SCHED("pumpMessages, iter through scheduledSends %p", *iter);
    if (find(sendsInFlight.begin(), sendsInFlight.end(), *iter) ==
        sendsInFlight.end()) {
      sendsInFlight.push_back(*iter);
      DEBUG_SCHED("%p calling propagate on first half", *iter);
      propagateBlkMsg(*iter);
      iter = scheduledSends.erase(iter);
    }
  }

  if (sendsInFlight.size() != 0 || scheduledSends.size() != 0) {
    //CkPrintf("Setting up idle condition; %d messages in flight\n", sendsInFlight.size());
    CcdCallOnCondition(CcdPROCESSOR_STILL_IDLE, pumpOnIdle, this);
  } else {
    DEBUG_SCHED("finished all sends");
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

void BlockScheduler::repositionBlock(StateList::iterator block) {
  StateList::iterator pos = block->pendingDependencies == 0 ?
    localBlocks.begin() : localBlocks.end();
  localBlocks.splice(pos, localBlocks, block);
}

template <typename K>
void BlockScheduler::updatePanel(map<K, Panel> &panels, K index) {
  typename map<K, Panel>::iterator iter = panels.find(index);
  CkAssert(iter != panels.end());
  Panel &panel = iter->second;

  panel.updatesLeftToPlan--;
  if(panel.updatesLeftToPlan == 0) {
    for (list<StateList::iterator>::iterator i = panel.dependents.begin();
         i != panel.dependents.end(); ++i) {
      (*i)->pendingDependencies--;
      repositionBlock(*i);
    }

    panels.erase(iter);
  }
}

template <typename K>
void BlockScheduler::addDependence(map<K, Panel> &panels, K index,
				   StateList::iterator block) {
  typename map<K, Panel>::iterator panel = panels.find(index);
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

  if (block.refs.size() == 1 && !block.isSending) {
    // First reference to this block, so ask for it
    DEBUG_SCHED("requesting getBlock from (%d, %d)", srcx, srcy);
    CkEntryOptions opts;
    luArr(src.first, src.second).
      getBlock(CkMyPe(), update->target->ix, update->target->iy,
               &(mgr->setPrio(GET_BLOCK, opts)));
  }

  if (block.m) {
    DEBUG_SCHED("already ARRIVED from (%d, %d)", srcx, srcy);
    update->tryDeliver(srcx, srcy, block.data);
  }
}

void BlockScheduler::deliverBlock(blkMsg *m) {
  DEBUG_SCHED("deliverBlock src (%d, %d)", m->src.x, m->src.y);

  CkAssert(wantedBlocks.find(make_pair(m->src)) != wantedBlocks.end());
  wantedBlock &block = wantedBlocks[make_pair(m->src)];
  block.m = m;
  block.data = m->data;

  for (list<Update*>::iterator update = block.refs.begin();
       update != block.refs.end(); ++update) {
    DEBUG_SCHED("tryDeliver of (%d,%d) to (%d,%d) @ %d", m->src.x, m->src.y,
		(*update)->target->ix, (*update)->target->iy, (*update)->t);
    (*update)->tryDeliver(m->src.x, m->src.y, m->data);
  }

  CkAssert(m->npes_receiver >= 1);
  CkAssert(CkMyPe() == m->pes[m->offset]);

  // This processor is no longer part of the set
  m->offset++;
  m->npes_receiver--;
  m->firstHalfSent = false;
  if (m->npes_receiver >= 1) {
    scheduleSend(m);
    block.isSending = true;
  }

  progress();
}

void BlockScheduler::propagateBlkMsg(blkMsg *m) {
  envelope *env = UsrToEnv(m);
  _SET_USED(env, 0);
  if (env->isPacked()) {
    unsigned char msgidx = env->getMsgIdx();
    if(_msgTable[msgidx]->unpack) {
      m = (blkMsg*)_msgTable[msgidx]->unpack(m);
      UsrToEnv(m)->setPacked(0);
    }
  }

  DEBUG_SCHED("(%d, %d): propagateBlkMsg npes = %p, firstHalfSent = %s",
              m->src.x, m->src.y, m->pes, m->firstHalfSent ? "true" : "false");

  if (!m->firstHalfSent) {
    LUBlk *block = luArr[m->src].ckLocal();
    if (block != NULL) {
      block->setupMsg(reverseSends);
      reverseSends = !reverseSends;
    }
    m->npes_sender = m->npes_receiver;
    m->npes_receiver = (m->npes_receiver+1) / 2;
  } else {
    int secondHalf = m->npes_sender - m->npes_receiver;
    CkAssert(secondHalf >= 0);
    m->offset += m->npes_receiver;
    m->npes_receiver = secondHalf;
  }

  DEBUG_SCHED("Delivering to processors, offset = %d, num = %d",
              m->offset, m->npes_receiver);
  for (int i = 0; i < m->npes_receiver; ++i)
    DEBUG_SCHED("\t%d", m->pes[i + m->offset]);

  if (m->npes_receiver >= 1) {
    takeRef(m);
    thisProxy[m->pes[m->offset]].deliverBlock(m);
  }

  traceMemoryUsage();
}

void BlockScheduler::dropRef(int srcx, int srcy, Update *update) {
  map<pair<int, int>, wantedBlock>::iterator input =
    wantedBlocks.find(make_pair(srcx, srcy));

  // Ignore local blocks used as inputs
  if (input == wantedBlocks.end())
    return;

  wantedBlock &block = input->second;
  block.refs.remove(update);
  if (block.refs.size() == 0 && !block.isSending) {
    delete block.m;
    wantedBlocks.erase(input);
  }
}

void BlockScheduler::runUpdate(list<Update>::iterator iter) {
  Update &update = *iter;
  CkAssert(update.ready());
  int tx = update.target->ix, ty = update.target->iy;
  update.triggered = true;
  pendingTriggered++;

  CkEntryOptions opts;
  int t = update.t;
  intptr_t update_ptr = (intptr_t)&update;
  luArr(tx, ty).processTrailingUpdate(t, update_ptr,
                                      &(mgr->setPrio(PROCESS_TRAILING_UPDATE,
                                                     opts, ty, tx, t)));
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
    //doneBlocks.erase(find(doneBlocks.begin(), doneBlocks.end(), *update.target));
  }

  pendingTriggered--;
  CkAssert(pendingTriggered >= 0);

  plannedUpdates.erase(find(plannedUpdates.begin(), plannedUpdates.end(), update));

  progress();
}

void BlockScheduler::updateUntriggered() {
  pendingTriggered--;
  if (pendingTriggered == 0) {
    pumpMessages();
  }
}

bool BlockScheduler::shouldExecute() {
  return numActive == 0;
}

void BlockScheduler::factorizationDone(CkIndex2D index) {
  DEBUG_SCHED("factorizationDone on (%d,%d)", index.x, index.y);

  if (index.x >= index.y)
    numActive--;

  map<pair<int, int>, list<Update*> >::iterator wanters =
    localWantedBlocks.find(make_pair(index.x, index.y));
  if (wanters != localWantedBlocks.end()) {
    list<Update*> &wantList = wanters->second;
    CkAssert(luArr(index).ckLocal());

    for (list<Update*>::iterator wanter = wantList.begin();
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
      for (list<ComputeU>::iterator computeU = pendingComputeU.begin();
           computeU != pendingComputeU.end(); ++computeU) {
        CkEntryOptions opts;
        luArr(computeU->x, computeU->y).
          processComputeU(0, &(mgr->setPrio(RECVL, opts, computeU->y)));
        pendingTriggered++;
        computeU = pendingComputeU.erase(computeU);
      }

      // Start trailing updates
      for (list<Update>::iterator update = plannedUpdates.begin();
	   update != plannedUpdates.end(); ++update) {
	if (update->ready()) {
          runUpdate(update);
          stateModified = true;
          break;
        }
      }
    }
  } while (stateModified);

  pumpMessages();

  inProgress = false;
}

