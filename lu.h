#include "manager.h"
#include "messages.h"
#include "scheduler.h"
#include "luConfig.h"
#include "lu.decl.h"
#include <ckmulticast.h>
#include <limits>
#include <map>
#include <algorithm>
#include <vector>
using std::min;

#if CHARMLU_DEBUG >= 1
    #define DEBUG_PRINT(FORMAT, ...) CkPrintf("(%d: [%d,%d]@%d) " FORMAT "\n", CkMyPe(), thisIndex.x, thisIndex.y, internalStep ,##__VA_ARGS__)
    #define DEBUG_SCHED(FORMAT, ...) CkPrintf("(%d S): " FORMAT "\n", CkMyPe() ,##__VA_ARGS__)
#else
    #define DEBUG_PRINT(...)
    #define DEBUG_SCHED(...)
#endif

#if CHARMLU_DEBUG >= 2
    #define DEBUG_PIVOT(...) CkPrintf(__VA_ARGS__)
    #define VERBOSE_PROGRESS(...) CkPrintf(__VA_ARGS__)
    #define VERY_VERBOSE_PIVOT_AGGLOM(...) CkPrintf(__VA_ARGS__)
    #define VERBOSE_VALIDATION(...) CkPrintf(__VA_ARGS__)
    #define VERBOSE_PIVOT_RECORDING
    #define VERBOSE_PIVOT_AGGLOM
#else
    #define DEBUG_PIVOT(...)
    #define VERBOSE_PROGRESS(...)
    #define VERY_VERBOSE_PIVOT_AGGLOM(...)
    #define VERBOSE_VALIDATION(...)
#endif

struct locval {
  locval(): val(0.0), loc(-1) {}
  locval(double _val, int _loc): val(_val), loc(_loc) {}
  double val;
  int loc;
};

/* readonly */
extern CProxy_Main mainProxy;
extern int traceTrailingUpdate;
extern int traceComputeU;
extern int traceComputeL;
extern int traceSolveLocalLU;
extern CkGroupID mcastMgrGID;

/// Global that holds the reducer type for locval
extern CkReduction::reducerType LocValReducer;
//extern CmiNodeLock lock;

static inline void takeRef(void *m) {
//    CmiLock(lock);
    CmiReference(UsrToEnv(m));
//    CmiUnlock(lock);
}
static inline void dropRef(void *m) {
//    CmiLock(lock);
    CmiFree(UsrToEnv(m));
//    CmiUnlock(lock);
}

class LUBlk: public CBase_LUBlk {
public:
  int internalStep;
  bool factored;
private:
  CProxy_BlockScheduler scheduler;
  BlockScheduler *localScheduler;
  int l_block, u_block;

  /// configuration settings
  LUConfig cfg;

  // The section of chares in the array on and below the current diagonal
  CProxySection_LUBlk belowLeft, belowRight;

  /// Variables used during factorization
  double *LU;

  int BLKSIZE, numBlks;
  blkMsg *L, *U;
  //BlockReadyMsg *mL, *mU;
  int activeCol, ind;

  LUMgr *mgr;

  /// Variables used only during solution
  double *bvec;
  //VALIDATION: variable to hold copy of untouched b vector (allocated during validation)
  double *b;
  //VALIDATION: variable to hold Ax
  double *Ax;

  double* storedVec;
  int diagRec;

  //VALIDATION: count of the number of point-to-point messages rcvd in a row
  int msgsRecvd;

  //VALIDATION: seed value used to regenerate A and b for validation
  int seed_A;
  int seed_b;

  // Timer for each block
  double startTime;

  //Variables for pivoting SDAG code

  int row1Index, row2Index, localRow1, localRow2,
    otherRowIndex, thisLocalRow, globalThisRow, globalOtherRow;
  bool remoteSwap;
  // Stores the local column max which is a candidate for that column's pivot element
  locval pivotCandidate;
  int pivotBlk;
  /// Tag for all msgs associated with a single batch of pivots
  int pivotBatchTag;
  /// A map to store and optimize the pivot operations
  std::map<int,int> pivotRecords;
  /// The number of rows factored since the last batch of pivots
  int numRowsSinceLastPivotSend ;
  /// The number of pending incoming remote pivots in a given batch
  int pendingIncomingPivots;
  /// The suggested pivot batch size
  int suggestedPivotBatchSize;
  /// How many blocks in the trailing submatrix have pulled our data and consumed it?
  int blockPulled;
  /// How many blocks in the trailing submatrix live on processors that have asked for it?
  int blocksAfter;
  /// Which PE's schedulers have requested this block?
  std::vector<CkArrayIndex1D> requestingPEs;

  /// The sub-diagonal chare array section that will participate in pivot selection
  /// @note: Only the diagonal chares will create and mcast along this section
  CProxySection_LUBlk pivotSection;
  /// All pivot sections members will save a cookie to their section
  CkSectionInfo pivotCookie;

  // The panel of blocks below the active diagonal chare
  CProxySection_LUBlk activePanel;

  /// The left-of-diagonal section of the chare array for pivoting
  CProxySection_LUBlk pivotLeftSection;

  /// The right-of-diagonal section of the chare array for pivoting
  CProxySection_LUBlk pivotRightSection;

  /// The optimal section for multicastRequestedBlock
  CProxySection_BlockScheduler panelAfter;
  std::set<int> panelAfterPEs;

  CProxySection_LUBlk rowBeforeDiag;
  CProxySection_LUBlk rowAfterDiag;
  CkSectionInfo rowBeforeCookie;
  CkSectionInfo rowAfterCookie;

  /// A pointer to the local branch of the multicast manager group that handles the pivot section comm
  CkMulticastMgr *mcastMgr;

  /// Pointer to a U msg indicating a pending L sub-block update. Used only in L chares
  UMsg *pendingUmsg;

  // For genBlock
  CrnStream blockStream, vecStream;

  LUBlk_SDAG_CODE

  public:
  LUBlk()
    : factored(false), storedVec(NULL), diagRec(0), msgsRecvd(0), blockPulled(0), blocksAfter(0)
  {
    __sdag_init();
#if defined(LU_TRACING)
    traceEnd();
#endif
  }

  void traceOn() {
#if defined(LU_TRACING)
    traceBegin();
#endif
  }

  //VALIDATION
  void startValidation();
  void recvXvec(int size, double* xvec);
  void sumBvec(int size, double* partial_b);
  void calcResiduals();
  void flushLogs();
  void finishInit();
  void initVec();
  void genBlock();
  void genVec(double *buf);
  void init(const LUConfig _cfg, CProxy_LUMgr _mgr, CProxy_BlockScheduler bs);
  void prepareForActivePanel(rednSetupMsg *msg);
  ~LUBlk();
  LUBlk(CkMigrateMessage* m) {}
  //added for migration
  void pup(PUP::er &p) {  }
  void computeU(blkMsg *givenLMsg);
  void computeL(blkMsg *givenUMsg);
  void updateMatrix(double *incomingL, double *incomingU);
  //broadcast the U downwards to the blocks in the same column
  void multicastRequestedBlock(PrioType prio);
  void multicastRecvU();
  void recvU(blkMsg *);
  //broadcast the L rightwards to the blocks in the same row
  void multicastRecvL();
  void sendBlocks(int);
  void getBlock(int pe);
  double *getBlock();
  int getActivePanelProgress();
  void processComputeU(int ignoredParam);
  void localSolve(double *xvec, double *preVec);
  void localForward(double *xvec);
  void localBackward(double *xvec);
  void beginForward(int size, double *preVec);
  void beginBackward(int size, double *preVec);
  void print();
  void print(const char* step);
private:
  // Copy received pivot data into its place in this block
  void applySwap(int row, int offset, double *data, double b);
  // Exchange local data
  void swapLocal(int row1, int row2, int offset=0);
  void doPivotLocal(int row1, int row2);
  /// Record the effect of a pivot operation in terms of actual row numbers
  void recordPivot(const int r1, const int r2);
  /** Is it time to send out the next batch of pivots
   *
   * @note: Any runtime adaptivity should be plugged here
   */
  bool shouldSendPivots();
  /// Periodically send out the agglomerated pivot operations
  void announceAgglomeratedPivots();
  /// Given a set of pivot ops, send out participating row chunks that you own
  void sendPendingPivots(const pivotSequencesMsg *msg);
  //internal functions for creating messages to encapsulate the priority
  inline blkMsg* createABlkMsg();
  locval findLocVal(int startRow, int col, locval first = locval());
  /// Update the sub-block of this L block starting at specified
  /// offset from the active column
  void updateLsubBlock(int activeCol, double* U, int offset=1, int startingRow=0);
  /// Compute the multipliers based on the pivot value in the
  /// received row of U and also find the candidate pivot in
  /// the immediate next column (after updating it simultaneously)
  locval computeMultipliersAndFindColMax(int col, double *U, int startingRow=0);
  inline int getIndex(int i, int j) {
    return i * BLKSIZE + j;
  }
};
