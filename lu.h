#include <charm++.h>
#include <vector>
#include <list>
#include <utility>

class BlockReadyMsg;

class BlockScheduler : public CBase_BlockScheduler {
public:
  BlockScheduler() : lock(CmiCreateLock()) { }

    void wantBlocks(CkIndex2D requester, BlockReadyMsg *mL, BlockReadyMsg *mU, int step);

  double* block(int block_id) { return blocks_ready[block_id].second; }
  void updateDone(CkIndex2D requester, int step) {}
  void blockArrived(int block_id) {}

private:
  CmiNodeLock lock;
  CProxy_LUBlk luArr;
  struct request {
    CkIndex2D requester;
    BlockReadyMsg *mL, *mU;
    int idxL, idxU;
    int step;
    request(CkIndex2D requester_, BlockReadyMsg *mL_, BlockReadyMsg *mU_, int step_)
      : requester(requester_), mL(mL_), mU(mU_), idxL(-1), idxU(-1), step(step_)
    { }
  };

  void processing(request r);
  void progress() {}

  std::list<request> reqs_waiting, reqs_pending, reqs_processing;
  /// Pointer to the data, and whether it's owned by a local LUBlk object
  std::vector<std::pair<bool, double *> > blocks_ready;
  std::vector<double *> blocks_free;
};
