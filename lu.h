#include <charm++.h>
#include <vector>
#include <list>
#include <utility>

class BlockReadyMsg;

class BlockScheduler : public CBase_BlockScheduler {
public:
    BlockScheduler(CProxy_LUBlk luArr_)
	: lock(CmiCreateLock()), luArr(luArr_)
	{ }

    void wantBlocks(CkIndex2D requester, BlockReadyMsg *mL, BlockReadyMsg *mU, int step);

  void updateDone(CkIndex2D requester, int step);
  void blockArrived(int x, int y);
  double *block(int x, int y);

private:
  CmiNodeLock lock;
  CProxy_LUBlk luArr;
  struct request {
    CkIndex2D requester;
    BlockReadyMsg *mL, *mU;
    int step;
    request(CkIndex2D requester_, BlockReadyMsg *mL_, BlockReadyMsg *mU_, int step_)
      : requester(requester_), mL(mL_), mU(mU_), step(step_)
    { }
  };

  struct block_state {
      double *data;
      int interested;

  block_state(double *d) : data(d), interested(0) {}
  };

  /// Drop a reference to the named block, and return whether this freed space 
  bool drop_block(int x, int y);

  void processing(request r);
  void progress();

  typedef std::list<request> req_list;
  req_list reqs_waiting, reqs_pending, reqs_processing;

  typedef std::map<std::pair<int, int>, block_state> block_map;
  block_map blocks_ready;
  std::vector<double *> blocks_free;
};
