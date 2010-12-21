/**
   @file A Charm++ implementation of LU 
   
   Messages of type blkMsg contain blocks of data that are sent down and to the 
   right as the computation progresses. 

   As the blkMsg messages arrive, they are always put into a buffer, and a call 
   to thisProxy(thisIndex.x,thisIndex.y).progress() is made if the computation 
   can proceed without any additional messages.

 */

#include "luConfig.h"

//#include <assert.h>
//#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <map>
#include <algorithm>
using std::min;
//#include <pthread.h>
#include <limits>

#if USE_CBLAS_H
extern "C" {
#include <cblas.h>
#include <clapack.h>
}

#elif USE_MKL_CBLAS_H
#include "mkl_cblas.h"
#include "mkl_lapack.h"

#elif USE_ACML_H
#include "acml.h"

#elif USE_ACCELERATE_BLAS
#include <Accelerate/Accelerate.h>

#elif USE_ESSL
#define _ESVCPTR
#include <complex>
#include <essl.h>


#else
#error "No BLAS Header files included!"
#endif

#ifndef ADAPT_SCHED_MEM
#error "Adaptive memory-aware scheduling is not enabled!"
#endif

#if USE_MEMALIGN
#include <malloc.h>
#endif

#include <comlib.h>
#include <controlPoints.h> // must come before user decl.h if they are using the pathInformationMsg
#include "lu.decl.h"
#include <trace-projections.h>
#include <ckmulticast.h>
#include <queueing.h> // for access to memory threshold setting

/* readonly: */
CProxy_Main mainProxy;
int traceTrailingUpdate;
int traceComputeU;
int traceComputeL;
int traceSolveLocalLU;
bool doPrioritize;
ComlibInstanceHandle multicastStats[4];
CProxy_locker lg;

#ifdef CHARMLU_DEBUG
    #define DEBUG_PRINT(...) CkPrintf(__VA_ARGS__)
    #define DEBUG_PIVOT(...) CkPrintf(__VA_ARGS__)
    #define VERY_VERBOSE_PIVOT_AGGLOM(...) CkPrintf(__VA_ARGS__)
    #define VERBOSE_PIVOT_RECORDING
    #define VERBOSE_PIVOT_AGGLOM
#else
    #define DEBUG_PRINT(...)
    #define DEBUG_PIVOT(...)
    #define VERY_VERBOSE_PIVOT_AGGLOM(...)
#endif

#include <cmath>

class locval {
    public:
        locval(): val(0.0), loc(-1) {}
        locval(double _val, int _loc): val(_val), loc(_loc) {}
        double val;
        int loc;
};

CkReductionMsg *maxLocVal(int nMsg, CkReductionMsg **msgs)
{
  CkAssert(nMsg > 0);

  locval *l = (locval*) msgs[0]->getData();
  for (int i = 1; i < nMsg; ++i) {
      locval *n = (locval *) msgs[i]->getData();
      if ( fabs(n->val) > fabs(l->val) )
          l = n;
  }

  return CkReductionMsg::buildNew(sizeof(locval), l);
}

/// Global that holds the reducer type for locval
CkReduction::reducerType LocValReducer;

/// Function that registers this reducer type on every processor
void registerLocValReducer()
{ LocValReducer = CkReduction::addReducer(maxLocVal); }

double infNorm(int size, double * array)
{
	  double maxval = fabs(array[0]);
	  for(int i = 1; i < size; i++) {
	      if(fabs(array[i]) > maxval)
	          maxval = fabs(array[i]);
	  }
	  return maxval;
}

//A class for randomly generating matrix elements' value
#define MAXINT (~(1<<31))
class MatGen {
private:
  //variables for generating the sequence of random numbers
  int curRnd;
  int rndA;
  int rndQ;
  int rndR;

public:
  MatGen(int seed) {
    curRnd = seed;
    rndA = 48271;
    rndQ = MAXINT/rndA;
    rndR = MAXINT%rndA;
  }

  int nextRndInt() {
    curRnd = rndA*(curRnd%rndQ) - rndR*(curRnd/rndQ);
    if (curRnd<0) curRnd += MAXINT;
    return curRnd;
  }

  //The range of the returned double random number is [-0.5, 0.5]
  double toRndDouble(int rndInt) {
    return (double)rndInt/MAXINT - 0.5;
  }

  double nextRndDouble() {
    return toRndDouble(nextRndInt());
  }

  void getNRndInts(int num, int *d) {
    for (int i=0; i<num; i++)
      d[i] = nextRndInt();
  }

  void getNRndDoubles(int num, double *d) {
    for (int i=0; i<num; i++)
      d[i] = nextRndDouble();
  }

  void skipNDoubles(int num) {
    for (int i=0; i<num; i++)
      nextRndDouble();
  }
};

struct blkMsg: public CMessage_blkMsg {
  // TODO: what is happening?
  char pad[16-((sizeof(envelope)+sizeof(int))%16)];
  double *data;

  void setMsgData(double *data_, int step, int BLKSIZE) {
    memcpy(data, data_, sizeof(double)*BLKSIZE*BLKSIZE);
    CkSetRefNum(this, step);
  }
};

class UMsg : public CMessage_UMsg, public CkMcastBaseMsg {
    public:
        UMsg(const int numElements, double *useg) {
            memcpy(data, useg, numElements * sizeof(double));
        }
        /// The post-diagonal portion of the active row
        double *data;
};

struct pivotMsg: public CMessage_pivotMsg, public CkMcastBaseMsg {
    /// @todo: active row is the same as row 1, right? Delete if not needed
  int activeRow, row1, row2;

  pivotMsg(int activeRow_, int row1_, int row2_) :
    activeRow(activeRow_), row1(row1_), row2(row2_) {
    CkSetRefNum(this, activeRow);
  }
};


class pivotSequencesMsg: public CMessage_pivotSequencesMsg, public CkMcastBaseMsg
{
    public:
        int numRowsProcessed;
        int numSequences;
        int *seqIndex;
        int *pivotSequence;

        pivotSequencesMsg(const int _firstRowProcessed, const int _numRowsProcessed)
        {
            numRowsProcessed = _numRowsProcessed;
            numSequences     = 0;
            CkSetRefNum(this, _firstRowProcessed);
        }
};



class pivotRowsMsg: public CMessage_pivotRowsMsg, public CkMcastBaseMsg
{
    public:
        int nRows;
        int blockSize;
        int *rowNum;
        double *rows;
        double *rhs;


        pivotRowsMsg(const int _blockSize, const int _refNum):
            nRows(0), blockSize(_blockSize)
        {
            CkSetRefNum(this, _refNum);
        }


        void copyRow(const int rNum, const double *row, const double b)
        {
            rowNum[nRows]    = rNum;
            rhs[nRows]       = b;
            memcpy(&rows[nRows*blockSize], row, blockSize*sizeof(double));
            nRows++;
        }
};



#include "manager.h"

class rednSetupMsg: public CkMcastBaseMsg, public CMessage_rednSetupMsg
{
    public:
        CkGroupID rednMgrGID;
        rednSetupMsg(CkGroupID _gid): rednMgrGID(_gid) {}
};

#include "mapping.h"

struct traceLU {
  int step, event;
  double startTime;
  traceLU(int internalStep, int eventType)
    : step(internalStep), event(eventType), startTime(CkWallTimer()) {
    traceUserSuppliedData(internalStep);
    traceMemoryUsage();
  }

  ~traceLU() {
    traceUserBracketEvent(event, startTime, CkWallTimer());
  }
};

CmiNodeLock lock;

struct locker : public CBase_locker {
    locker() { lock = CmiCreateLock(); }
};

static inline void takeRef(void *m) {
    CmiLock(lock);
    CmiReference(UsrToEnv(m));
    CmiUnlock(lock);
}
static inline void dropRef(void *m) {
    CmiLock(lock);
    CmiFree(UsrToEnv(m));
    CmiUnlock(lock);
}

class Main : public CBase_Main {
  LUConfig luCfg;
  double startTime;
  int iteration;
  int numIterations;

  int gMatSize;
  int whichMulticastStrategy;
  bool solved, LUcomplete, workStarted;
  bool sentVectorData;

  CProxy_LUBlk luArrProxy;

public:
    Main(CkArgMsg* m) : iteration(0), numIterations(1), solved(false), LUcomplete(false), workStarted(false), sentVectorData(false) {

    if (m->argc<4) {
      CkPrintf("Usage: %s <matrix size> <block size> <mem threshold> [<pivot batch size> <mapping scheme> [<peTileRows> <peTileCols>] ]\n", m->argv[0]);
      CkExit();
    }

    /// Parse the command line and accept user input
    gMatSize             = atoi(m->argv[1]);
    luCfg.blockSize      = atoi(m->argv[2]);
    luCfg.memThreshold   = atoi(m->argv[3]);

    if (m->argc >= 5)
        luCfg.pivotBatchSize = atoi( m->argv[4] );
    else
        luCfg.pivotBatchSize = luCfg.blockSize / 4;

    if (m->argc >= 6)
    {
        luCfg.mappingScheme = atoi( m->argv[5] );
        if (luCfg.mappingScheme == 3)
        {
            if (m->argc < 8) {
                std::pair<int,int> tileDims = computePETileDimensions();
                luCfg.peTileRows = tileDims.first;
                luCfg.peTileCols = tileDims.second;
            }
            else {
                luCfg.peTileRows = atoi( m->argv[6] );
                luCfg.peTileCols = atoi( m->argv[7] );
            }
            int peTileSize = luCfg.peTileRows * luCfg.peTileCols;
            if ( peTileSize > CkNumPes() )
                CkAbort("The PE tile dimensions are too big for the num of PEs available!");
            if ( peTileSize < CkNumPes() )
                CkPrintf("WARNING: Configured to use a PE tile size(%dx%d) (for 2D tile mapping) that does not use all the PEs(%d)\n",
                         luCfg.peTileRows, luCfg.peTileCols, CkNumPes() );
        }
    }

    numIterations = 1;
    CkPrintf("CLI: numIterations=%d\n", numIterations);

    if (gMatSize%luCfg.blockSize!=0) {
      CkPrintf("The matrix size %d should be a multiple of block size %d!\n", gMatSize, luCfg.blockSize);
      CkExit();
    }

    luCfg.numBlocks = gMatSize / luCfg.blockSize;

    CkPrintf("Running LU on %d processors (%d nodes): "
             "\n\tMatrix size: %d X %d "
             "\n\tBlock size: %d X %d "
             "\n\tPivot batch size: %d"
             "\n\tMem Threshold (MB): %d"
             "\n\tMapping Scheme: %d (%s)\n",
             CkNumPes(), CmiNumNodes(),
             gMatSize, gMatSize,
             luCfg.blockSize, luCfg.blockSize,
             luCfg.pivotBatchSize,
             luCfg.memThreshold,
             luCfg.mappingScheme,
             luCfg.mappingScheme == 1 ? "Balanced Snake" : (luCfg.mappingScheme==2? "Block Cylic": "2D Tiling")
             );
    if (luCfg.mappingScheme == 3)
        CkPrintf("\tMapping PE tile size: %d x %d\n", luCfg.peTileRows, luCfg.peTileCols);



    mainProxy = thisProxy;
    doPrioritize = false;
      
    DEBUG_PRINT("Registering user events\n");
    traceTrailingUpdate = traceRegisterUserEvent("Trailing Update");
    traceComputeU = traceRegisterUserEvent("Compute U");
    traceComputeL = traceRegisterUserEvent("Compute L");
    traceSolveLocalLU = traceRegisterUserEvent("Solve local LU");
    
    traceRegisterUserEvent("Local Multicast Deliveries", 10000);    
    traceRegisterUserEvent("Remote Multicast Forwarding - preparing", 10001);
    traceRegisterUserEvent("Remote Multicast Forwarding - sends", 10002);

    multicastStats[0] = ComlibRegister(new OneTimeRingMulticastStrategy() ); 
    multicastStats[1] = ComlibRegister(new OneTimeNodeTreeMulticastStrategy(2) ); 
    multicastStats[2] = ComlibRegister(new OneTimeNodeTreeMulticastStrategy(3) ); 
    multicastStats[3] = ComlibRegister(new OneTimeNodeTreeMulticastStrategy(4) ); 

    ControlPoint::EffectIncrease::GrainSize("block_size");
    ControlPoint::EffectDecrease::Concurrency("block_size");
    
    //	  ControlPoint::EffectIncrease::Concurrency("mapping");
    ControlPoint::EffectIncrease::NumMessages("mapping");
    ControlPoint::EffectIncrease::MessageOverhead("mapping");

    ControlPoint::EffectIncrease::MemoryConsumption("memory_threshold");

    lg = CProxy_locker::ckNew();

    thisProxy.iterationCompleted();
  }

  void finishInit() {
    if (!sentVectorData) {
      luArrProxy.initVec();
      sentVectorData = true;
    } else {
      luArrProxy.flushLogs();
    }
  }

  void continueIter() {
    CkCallback *cb = new CkCallback(CkIndex_Main::arrayIsCreated(NULL), thisProxy);
    CkStartQD(*cb); // required currently for use with new Comlib
  }

  void arrayIsCreated(CkReductionMsg * m) {
    startTime = CmiWallTimer();
    luArrProxy.factor();
  }

  void iterationCompleted() {
    CkPrintf("called iterationCompleted()\n");

    if (workStarted) {
      LUcomplete = true;
    }

    if (solved && LUcomplete) {
      double endTime = CmiWallTimer();
      double duration = endTime - startTime;
      registerControlPointTiming(duration);
      
      CkPrintf("Iteration %d time: %fs\n", iteration, duration);
      
      outputStats();

      //VALIDATION: Instead of ending program, perform validation
      //terminateProg();
      luArrProxy.startValidation();
    } else if (!solved && LUcomplete) {
      luArrProxy.print();
      luArrProxy(0,0).forwardSolve();
      solved = true;
      iteration++;
    } else {
      // Only advance phases after a few factorizations have been performed
      // Prior to the first phase of actual work, iteration=1
      if( 1 || iteration % 2 == 1 || iteration==1){
#if 0
	gotoNextPhase();
      
	whichMulticastStrategy = controlPoint("multicast_strategy", 2, 2);
	BLKSIZE = 1 << 2; //1 << controlPoint("block_size", 10,10);
	mapping = controlPoint("mapping", 1, 1);
	memThreshold = 200 + controlPoint("memory_threshold", 0, 20) * 100;
      
	// CkPrintf("%d %d %d\n",  (int)BLKSIZE, (int)mapping, (int)whichMulticastStrategy);
	// fflush(stdout);
	// fflush(stderr);
      
	numBlks = gMatSize/BLKSIZE;
#endif
      }
    
    
      char note[200];
      sprintf(note, "*** New iteration: block size = %d, mapping = %d %s, multicast = %d, memthreshold = %d MB", 
	      luCfg.blockSize, luCfg.mappingScheme, luCfg.mappingScheme == 1 ? "Balanced Snake" : (luCfg.mappingScheme==2? "Block Cylic": "2D Tiling"), whichMulticastStrategy, luCfg.memThreshold);
      traceUserSuppliedNote(note);
      CkPrintf("%s\n", note);
      fflush(stdout);
    
      CkArrayOptions opts(luCfg.numBlocks, luCfg.numBlocks);
      opts.setAnytimeMigration(false)
	  .setStaticInsertion(true);
      switch (luCfg.mappingScheme) {
      case 0:
	opts.setMap(CProxy_BlockCyclicMap::ckNew());
	break;
      case 1:
        opts.setMap(CProxy_LUBalancedSnakeMap::ckNew(luCfg.numBlocks, luCfg.blockSize));
	break;
      case 2:
	  opts.setMap(CProxy_RealBlockCyclicMap::ckNew(1, luCfg.numBlocks));
	  break;
      case 3:
          opts.setMap( CProxy_PE2DTilingMap::ckNew(luCfg.peTileRows, luCfg.peTileCols) );
          break;
      }
      CProxy_LUMgr mgr = CProxy_PrioLU::ckNew(luCfg.blockSize, gMatSize);

      luArrProxy = CProxy_LUBlk::ckNew(opts);

      workStarted = true;
    
      luArrProxy.startup(luCfg, 0, mgr);
    }
  }


  std::pair<int,int> computePETileDimensions()
  {
      // Identify two factors that can be used as the tile dimensions for the PE tile
      int factor1 = sqrt(CkNumPes());
      while ( (CkNumPes() % factor1 != 0) && (factor1 > 0) )
          factor1--;
      if (factor1 == 0)
          CkAbort("Couldn't identify a factor of numPEs to represent the PEs as a 2D tile");

      int factor2 = CkNumPes() / factor1;

      // Set the tile dimensions
      int numPErows = (factor1 >= factor2) ? factor1 : factor2;
      int numPEcols = CkNumPes() / numPErows;

      if (numPErows * numPEcols != CkNumPes())
          CkAbort("The identified tile dimensions dont match the number of PEs!!");

      return std::make_pair(numPErows, numPEcols);
  }

  
  void outputStats() {
    double endTime = CmiWallTimer();
    double duration = endTime-startTime;

    double n = gMatSize;

    long long flopCount = 0;	 // floating point ops
    for (int i=1;i<=gMatSize;i++) {
      for (int j=1+i; j<=gMatSize; j++) {
	flopCount += (1+2*gMatSize-2*i);
      }
    }

    double flops = ((double)flopCount)	/ duration; // floating point ops per second
    double gflops = flops / 1000000000.0; // Giga fp ops per second
    std::cout << "RESULT procs: \t" << CkNumPes() << "\tblock size:\t" << luCfg.blockSize << "\tGFlops:\t" << gflops << "\tTime(s):\t" << duration << std::endl;

    double HPL_flop_count =  (2.0/3.0*n*n*n+3.0/2.0*n*n)/duration ;
    double HPL_gflops =	 HPL_flop_count / 1000000000.0; // Giga fp ops per second
    std::cout << "HPL flop count gives \t" << HPL_gflops << "\tGFlops" << std::endl;

    
    double gflops_per_core = HPL_gflops / (double)CkNumPes();

    struct {
      const char *machine;
      double gflops_per_core;
    } peaks[] = {{"order.cs.uiuc.edu", 7.4585},
		 {"abe.ncsa.uiuc.edu", 9.332},
		 {"Kraken", 10.4},
		 {"BG/P", 3.4}};

    for (int i = 0; i < sizeof(peaks)/sizeof(peaks[0]); ++i) {
      double fractionOfPeak = gflops_per_core / peaks[i].gflops_per_core;
      std::cout << "If ran on " << peaks[i].machine << ", I think you got \t"
		<< 100.0*fractionOfPeak << "% of peak" << std::endl;
    }
  }

  void calcScaledResidual(CkReductionMsg *msg) {
	int reducedArrSize=msg->getSize() / sizeof(double);
	double *maxvals=(double *) msg->getData();

	double n = luCfg.blockSize * luCfg.numBlocks;

	double r = maxvals[3]/((maxvals[0]*maxvals[2]+maxvals[1])*n*std::numeric_limits<double>::epsilon());

	DEBUG_PRINT("|A|inf = %e\n|b|inf = %e\n|x|inf = %e\n|Ax-b|inf = %e\n",maxvals[0],maxvals[1],maxvals[2],maxvals[3]);

	CkPrintf("epsilon = %e\nresidual = %f\n",std::numeric_limits<double>::epsilon(),r);
	if(r>16)
		CkPrintf("=== WARNING: Scaled residual is greater than 16 - OUT OF SPEC ===\n");

	delete msg;

    CkCallback cb(CkIndex_Main::done(NULL),thisProxy); 
    traceCriticalPathBack(cb);
  }

  void done(pathInformationMsg *m){
    // CkPrintf("Main::done() After critical path has been determined\n");
    //	  m->printme();
    gotoNextPhase(); // Make sure we get timings for the phase that just finished.
    CkExit();
  }

};




class LUBlk: public CBase_LUBlk {

    /// configuration settings
    LUConfig cfg;

  // The section of chares in the array on and below the current diagonal
  CProxySection_LUBlk belowLeft, belowRight;

  /// Variables used during factorization
  double **LU;

  int BLKSIZE, numBlks;
  blkMsg *L, *U;
  int internalStep, activeCol, ind;

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

  CProxySection_LUBlk rowBeforeDiag;
  CProxySection_LUBlk rowAfterDiag;
  CkSectionInfo rowBeforeCookie;
  CkSectionInfo rowAfterCookie;

  /// A pointer to the local branch of the multicast manager group that handles the pivot section comm
  CkMulticastMgr *mcastMgr;

  /// Pointer to a U msg indicating a pending L sub-block update. Used only in L chares
  UMsg *pendingUmsg;

  LUBlk_SDAG_CODE

  // The following declarations are used for optimization and
  // analysis. They are not essential to the algorithm.
  int whichMulticastStrategy;

public:
  LUBlk() : storedVec(NULL), diagRec(0), msgsRecvd(0) {
      __sdag_init();
    whichMulticastStrategy = 0;

    //CkAssert(BLKSIZE>0); // If this fails, readonly variables aren't
			 // propagated soon enough. I'm assuming they
			 // are safe to use here.




    /*internalStep = 0;	 
     
    traceUserSuppliedData(-1);	
    traceMemoryUsage();	 
     
    MatGen rnd(thisIndex.x * numBlks + thisIndex.y);	  
    for (int i=0; i<BLKSIZE*BLKSIZE; i++) {  
      LU[i] = rnd.toRndDouble(rnd.nextRndInt());  
    } 


    testdgemm();*/

  }

  //VALIDATION
  void startValidation() {
	  // Starting state:
	  // solution sub-vector x is in variable bvec on the diagonals
	  // variable b has the original b vector on the diagonals

	  //Diagonals regenerate b and distribute x across entire column
	  if(thisIndex.x == thisIndex.y)
	  {

		  CProxySection_LUBlk col =
				  CProxySection_LUBlk::ckNew(thisArrayID, 0, numBlks-1,
						  1, thisIndex.y, thisIndex.y, 1);

		  col.recvXvec(BLKSIZE, bvec);
	  }
  }

  //VALIDATION
  void recvXvec(int size, double* xvec) {
      //Regenerate A and place into already allocated LU
      genBlock();


	  double *partial_b = new double[BLKSIZE];

	  //Perform local dgemv
#if USE_ESSL
    dgemv("T", BLKSIZE, BLKSIZE, 1.0, LU[0], BLKSIZE, xvec, 1, 0.0, partial_b, 1);
#else
    cblas_dgemv( CblasRowMajor, CblasNoTrans,
              BLKSIZE, BLKSIZE, 1.0, LU[0],
              BLKSIZE, xvec, 1, 0.0, partial_b, 1);
#endif

    //sum-reduction of result across row with diagonal element as target
    thisProxy(thisIndex.x,thisIndex.x).sumBvec(BLKSIZE,partial_b);

    //if you are not the diagonal, find your max A value and contribute
    if(thisIndex.x != thisIndex.y) {
    	//find local max of A
        double A_max = infNorm(BLKSIZE * BLKSIZE, LU[0]);
    	DEBUG_PRINT("[%d,%d] A_max  = %e\n",thisIndex.x,thisIndex.y,A_max);

    	double maxvals[4];
    	maxvals[0] = A_max;
    	maxvals[1] = -1;
    	maxvals[2] = -1;
    	maxvals[3] = -1;

    	contribute(sizeof(maxvals), &maxvals, CkReduction::max_double,CkCallback(CkIndex_Main::calcScaledResidual(NULL),mainProxy));
    }
  }

  //VALIDATION
  void sumBvec(int size, double* partial_b) {

	  //Clear bvec before first message processed for sum-reduction
	  if(msgsRecvd == 0) {
	      Ax = new double[BLKSIZE];
	      memset(Ax, 0, BLKSIZE*sizeof(double));
	  }

	  //Sum up messages
	  if(++msgsRecvd <= numBlks) {
		  for (int i = 0; i < size; i++) {
			  Ax[i] += partial_b[i];
		  }
	  }

	  //if all messages recieved, calculate the residual
	  if (msgsRecvd == numBlks)
		  calcResiduals();
  }

  //VALIDATION
  void calcResiduals() {
      b = new double[BLKSIZE];
      genVec(b);
	  double *residuals = new double[BLKSIZE];

	  //diagonal elements that received sum-reduction perform b - A*x
	  for (int i = 0; i < BLKSIZE; i++) {
		  residuals[i] = b[i] - Ax[i];
//		  if(fabs(residuals[i]) > 1e-14 || std::isnan(residuals[i]) || std::isinf(residuals[i]))
//			  CkPrintf("WARNING: Large Residual for x[%d]: %f - %f = %e\n", thisIndex.x*BLKSIZE+i, b[i], bvec[i], residuals[i]);
	  }

	  //find local max values
	  double A_max = infNorm(BLKSIZE * BLKSIZE, LU[0]);
	  double b_max = infNorm(BLKSIZE, b);
	  double x_max = infNorm(BLKSIZE, bvec);
	  double res_max = infNorm(BLKSIZE, residuals);
	  DEBUG_PRINT("[%d,%d] A_max  = %e\n",thisIndex.x,thisIndex.y,A_max);
	  DEBUG_PRINT("[%d,%d] b_max  = %e\n",thisIndex.x,thisIndex.y,b_max);
	  DEBUG_PRINT("[%d,%d] x_max  = %e\n",thisIndex.x,thisIndex.y,x_max);
	  DEBUG_PRINT("[%d,%d] res_max  = %e\n",thisIndex.x,thisIndex.y,res_max);

	  double maxvals[4];
	  maxvals[0] = A_max;
	  maxvals[1] = b_max;
	  maxvals[2] = x_max;
	  maxvals[3] = res_max;

	  contribute(sizeof(maxvals), &maxvals, CkReduction::max_double,CkCallback(CkIndex_Main::calcScaledResidual(NULL),mainProxy));
  }

  void flushLogs() {
    flushTraceLog();
    contribute(CkCallback(CkIndex_Main::continueIter(), mainProxy));	
  }

  void testdgemm(){
#if 0
    unsigned long blocksize = 32 << thisIndex.x;
    
    if(thisIndex.y == 0 && thisIndex.x == 1){	  
      
#if USE_MEMALIGN
      double *m1 = (double*)memalign(128, blocksize*blocksize*sizeof(double) );
      double *m2 = (double*)memalign(128, blocksize*blocksize*sizeof(double) ); 
      double *m3 = (double*)memalign(128, blocksize*blocksize*sizeof(double) );
      if(m1 == NULL || m2 == NULL || m3 == NULL)
	return;
#else
      double *m1 = new double[blocksize*blocksize]; 
      double *m2 = new double[blocksize*blocksize]; 
      double *m3 = new double[blocksize*blocksize]; 
#endif

      MatGen rnd(0); 

      rnd.getNRndDoubles(blocksize * blocksize, m1);
      rnd.getNRndDoubles(blocksize * blocksize, m2);
      rnd.getNRndDoubles(blocksize * blocksize, m3);

      double startTest = CmiWallTimer(); 
      
#if USE_ESSL
      dgemm( "N", "N",
	     blocksize, blocksize, blocksize,
	     -1.0, m1,
	     blocksize, m2, blocksize,
	     1.0, m3, blocksize);
#else
      cblas_dgemm( CblasRowMajor, 
		   CblasNoTrans, CblasNoTrans, 
		   blocksize, blocksize, blocksize, 
		   -1.0, m1, 
		   blocksize, m2, blocksize, 
		   1.0, m3, blocksize); 
#endif	   
      
      double endTest = CmiWallTimer(); 
      double duration = endTest-startTest; 
      
      CkPrintf("The dgemm %d x %d call takes %g seconds\n", blocksize, blocksize, duration); 
      double flopcount = (double)blocksize * (double)blocksize * (double)blocksize * 2.0; 
      double gflopcount = flopcount / 1000000000.0; 
      double gflopPerSec = gflopcount / duration; 
      CkPrintf("The dgemm\t%d\tx %d achieves\t%g\tGFlop/sec\n", blocksize, blocksize, gflopPerSec); 

      {

#if USE_ESSL
      dgemm( "T", "T",
	     blocksize, blocksize, blocksize,
	     -1.0, m1,
	     blocksize, m2, blocksize,
	     1.0, m3, blocksize);
#else
      cblas_dgemm( CblasRowMajor, 
		   CblasTrans, CblasTrans, 
		   blocksize, blocksize, blocksize, 
		   -1.0, m1, 
		   blocksize, m2, blocksize, 
		   1.0, m3, blocksize); 
#endif	   
      
      double endTest = CmiWallTimer(); 
      double duration = endTest-startTest; 
      
      CkPrintf("The dgemm %d x %d Transpose call takes %g seconds\n", blocksize, blocksize, duration); 
      double flopcount = (double)blocksize * (double)blocksize * (double)blocksize * 2.0; 
      double gflopcount = flopcount / 1000000000.0; 
      double gflopPerSec = gflopcount / duration; 
      CkPrintf("The dgemm\t%d\tx %d Transpose achieves\t%g\tGFlop/sec\n", blocksize, blocksize, gflopPerSec); 
      
      }

      delete[] m1;
      delete[] m2;
      delete[] m3;
    } 
#endif
  }

  void initVec() {
    bvec = new double[BLKSIZE];

    seed_b = 9934835;
    genVec(bvec);

#if defined(PRINT_VECTORS)
    for (int i = 0; i < BLKSIZE; i++) {
      CkPrintf("memcpy bvec[%d] = %f\n", i, bvec[i]);
    }
#endif

    contribute(CkCallback(CkIndex_Main::finishInit(), mainProxy));
  }

    CrnStream blockStream, vecStream;

  void genBlock()
    {
	CrnStream stream;
	memcpy(&stream, &blockStream, sizeof(CrnStream));

	for (double *d = LU[0]; d < LU[0] + BLKSIZE*BLKSIZE; ++d)
	    *d = CrnDouble(&stream);
    }

  void genVec(double *buf)
    {
      MatGen rnd(seed_b);

      // Skip the blocks before this one
      rnd.skipNDoubles(thisIndex.x * BLKSIZE);
      rnd.getNRndDoubles(BLKSIZE, buf);
    }

  void init(const LUConfig _cfg, int _whichMulticastStrategy, CProxy_LUMgr _mgr) {
    cfg = _cfg;
    whichMulticastStrategy = _whichMulticastStrategy;
    BLKSIZE = cfg.blockSize;
    numBlks = cfg.numBlocks;
    mgr = _mgr.ckLocalBranch();
    suggestedPivotBatchSize = cfg.pivotBatchSize;
    
    // Set the schedulers memory usage threshold to the one based upon a control point
    schedAdaptMemThresholdMB = cfg.memThreshold;

    CkAssert(BLKSIZE>0); // If this fails, readonly variables aren't
			 // propagated soon enough. I'm assuming they
			 // are safe to use here.

#if USE_MEMALIGN
    LU    = (double**) memalign(128, BLKSIZE*sizeof(double*) );
    CkAssert(LU != NULL);
    LU[0] = (double*)  memalign(128, BLKSIZE*BLKSIZE*sizeof(double) );
    CkAssert(LU[0] != NULL);
#else
    LU    = new double* [BLKSIZE];
    LU[0] = new double  [BLKSIZE*BLKSIZE];
#endif

    // Initialize the pointers to the rows of the matrix block
    for (int i=1; i<BLKSIZE; i++)
        LU[i] = LU[0] + i*BLKSIZE;

    internalStep = 0;  
     
    traceUserSuppliedData(-1);	
    traceMemoryUsage();	 
     
    //VALIDATION: saved seed value to use for validation
    seed_A = 2998389;
    CrnInitStream(&blockStream, seed_A + thisIndex.x*numBlks + thisIndex.y, 0);
    genBlock();

#if 0
    double b = thisIndex.x * BLKSIZE + 1, c = thisIndex.y * BLKSIZE + 1;
    for (int i = 0; i<BLKSIZE*BLKSIZE; i++) {
      if (i % BLKSIZE == 0 && thisIndex.y == 0) {
	LU[i] = b;
	b += 1.0;
	c += 1.0;
      } else if (i < BLKSIZE && thisIndex.x == 0) {
	LU[i] = c;
	c += 1.0;
      } else 
	LU[i] = 0.0;
    }

    if (thisIndex.x == thisIndex.y) {
      b = thisIndex.x * BLKSIZE + 1;

      for (int i = 0; i < BLKSIZE*BLKSIZE; i+=BLKSIZE+1) {
	LU[i] = b;
	b += 1.0;
      }
    }
#endif

    this->print("input-generated-LU");

    testdgemm();

    /// Chares on the array diagonal will now create pivot sections that they will talk to
    if (thisIndex.x == thisIndex.y)
    {
        // Create the pivot section
        pivotSection = CProxySection_LUBlk::ckNew(thisArrayID, thisIndex.x,numBlks-1,1,thisIndex.y,thisIndex.y,1);
        activePanel = CProxySection_LUBlk::ckNew(thisArrayID, thisIndex.x+1,numBlks-1,1,thisIndex.y,thisIndex.y,1);
        pivotLeftSection = CProxySection_LUBlk::ckNew(thisArrayID, thisIndex.x, numBlks-1, 1, 0, thisIndex.y-1, 1);
        pivotRightSection = CProxySection_LUBlk::ckNew(thisArrayID, thisIndex.x, numBlks-1, 1, thisIndex.y+1, numBlks-1, 1);
        rowBeforeDiag = CProxySection_LUBlk::ckNew(thisArrayID, thisIndex.x,thisIndex.x,1,0,thisIndex.y-1,1);
        rowAfterDiag = CProxySection_LUBlk::ckNew(thisArrayID, thisIndex.x,thisIndex.x,1,thisIndex.y+1,numBlks-1,1);
        // Create a multicast manager group
        CkGroupID mcastMgrGID = CProxy_CkMulticastMgr::ckNew();
        CkMulticastMgr *mcastMgr = CProxy_CkMulticastMgr(mcastMgrGID).ckLocalBranch();
        // Delegate pivot section to the manager
        pivotSection.ckSectionDelegate(mcastMgr);
        activePanel.ckSectionDelegate(mcastMgr);
        pivotLeftSection.ckSectionDelegate(mcastMgr);
        pivotRightSection.ckSectionDelegate(mcastMgr);
        rowBeforeDiag.ckSectionDelegate(mcastMgr);
        rowAfterDiag.ckSectionDelegate(mcastMgr);
        // Set the reduction client for this pivot section
        mcastMgr->setReductionClient( pivotSection, new CkCallback( CkIndex_LUBlk::colMax(0), thisProxy(thisIndex.y, thisIndex.y) ) );

        // Invoke a dummy mcast so that all the section members know which section to reduce along
        rednSetupMsg *pivotMsg = new rednSetupMsg(mcastMgrGID);
        rednSetupMsg *activePanelMsg = new rednSetupMsg(mcastMgrGID);
        rednSetupMsg *pivotLeftMsg = new rednSetupMsg(mcastMgrGID);
        rednSetupMsg *pivotRightMsg = new rednSetupMsg(mcastMgrGID);
        rednSetupMsg *rowBeforeMsg = new rednSetupMsg(mcastMgrGID);
        rednSetupMsg *rowAfterMsg = new rednSetupMsg(mcastMgrGID);

        pivotSection.prepareForPivotRedn(pivotMsg);
        activePanel.prepareForActivePanel(activePanelMsg);
        pivotLeftSection.prepareForPivotLR(pivotLeftMsg);
        pivotRightSection.prepareForPivotLR(pivotRightMsg);
        rowBeforeDiag.prepareForRowBeforeDiag(rowBeforeMsg);
        rowAfterDiag.prepareForRowAfterDiag(rowAfterMsg);

        if (thisIndex.x == 0) {
          thisProxy.multicastRedns(0);
        }
    }

    // All chares except members of pivot sections are done with init
  }

  ~LUBlk() {
    //CkPrintf("freeing LuBlk\n");
#if USE_MEMALIGN
    free(LU[0]);
    free(LU);
#else
    delete [] LU[0];
    delete [] LU;
#endif
    LU = NULL;
  }

  LUBlk(CkMigrateMessage* m) {}

  //added for migration
  void pup(PUP::er &p) {  }
 
//   // Called on element 0,0
//   void begin() {
//     CkAssert(thisIndex.x == 0 && thisIndex.y == 0);
//     traceUserSuppliedData(internalStep);
//     traceMemoryUsage();
//     thisProxy(0,0).processLocalLU(0);
//   }

  /* Computation functions that should be called localy related with each
   * state of the block
   */

  void computeU(blkMsg *givenLMsg) {
    traceLU t(internalStep, traceComputeU);
    double *givenL = givenLMsg->data;

#if 0
    if( ((unsigned long)givenL) % 16 != 0){
      CkPrintf("givenL mod 16=%d\n", (int)(((unsigned long)givenL) % 16 ));
      CkPrintf("sizeof(envelope)=%d\n", sizeof(envelope));
      CkPrintf("sizeof(int) = %d\n", sizeof(int));
      //      "16-((sizeof(envelope)+sizeof(int))%16";
    }

    CkAssert( ((unsigned long)givenL) % 16 == 0);
#endif

    DEBUG_PRINT("elem[%d,%d]::computeU called at step %d\n", thisIndex.x, thisIndex.y, internalStep);

    //processing row by row (forward substitution)
    //the 1st row of U is not changed

    //solve following rows based on previously solved rows
    //row indicates the row of U that is just solved

#if USE_ESSL
    // givenL is implicitly transposed by telling dtrsm that it is a
    // right, upper matrix. Since this also switches the order of
    // multiplication, the transpose is output to LU.
    dtrsm("R", "U", "N", "U", BLKSIZE, BLKSIZE, 1.0, givenL, BLKSIZE, LU[0], BLKSIZE);
#else
    cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, BLKSIZE, BLKSIZE, 1.0, givenL, BLKSIZE, LU[0], BLKSIZE);
#endif
  }

  void computeL(blkMsg *givenUMsg) {
    traceLU t(internalStep, traceComputeL);
    double *givenU = givenUMsg->data;
    //	  CkAssert( ((unsigned long)givenU) % 16 == 0);

    DEBUG_PRINT("elem[%d,%d]::computeL called at step %d\n", thisIndex.x, thisIndex.y, internalStep);

#if USE_ESSL
    dtrsm("R", "U", "N", "N", BLKSIZE, BLKSIZE, 1.0, givenU, BLKSIZE, LU[0], BLKSIZE);
#else
    cblas_dtrsm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, BLKSIZE, BLKSIZE, 1.0, givenU, BLKSIZE, LU[0], BLKSIZE);
#endif
  }

  void updateMatrix(blkMsg *givenLMsg, blkMsg *givenUMsg) {
    traceLU t(internalStep, traceTrailingUpdate);
    DEBUG_PRINT("elem[%d,%d] is updating its value at step %d\n", thisIndex.x, thisIndex.y, internalStep);

    double *incomingL = givenLMsg->data;
    double *incomingU = givenUMsg->data;

#if USE_ESSL
    // By switching the order of incomingU and incomingL the transpose
    // is applied implicitly: C' = B*A
    dgemm( "N", "N",
	   BLKSIZE, BLKSIZE, BLKSIZE,
	   -1.0, incomingU,
	   BLKSIZE, incomingL, BLKSIZE,
	   1.0, LU[0], BLKSIZE);
#else
    cblas_dgemm( CblasRowMajor,
		 CblasNoTrans, CblasNoTrans,
		 BLKSIZE, BLKSIZE, BLKSIZE,
		 -1.0, incomingL,
		 BLKSIZE, incomingU, BLKSIZE,
		 1.0, LU[0], BLKSIZE);
#endif
  }

  //broadcast the U downwards to the blocks in the same column
  inline void multicastRecvU() {
    traceUserSuppliedData(internalStep);
    traceMemoryUsage();
    
    DEBUG_PRINT("[PE %d] elem %d,%d Multicast to part of column %d from step %d\n", CkMyPe(), thisIndex.x, thisIndex.y, thisIndex.y, internalStep);
    
    CProxySection_LUBlk oneCol = CProxySection_LUBlk::ckNew(thisArrayID, thisIndex.x+1, numBlks-1, 1, thisIndex.y, thisIndex.y, 1);

    if (whichMulticastStrategy > -1)
      ComlibAssociateProxy(multicastStats[whichMulticastStrategy], oneCol);
    
    blkMsg *givenU = createABlkMsg();
    *(int*)CkPriorityPtr(givenU) = -1;
    oneCol.recvU(givenU);

//     for(int i=thisIndex.x+1; i<numBlks; i++){
//	 blkMsg *givenU = createABlkMsg();
//	 DEBUG_PRINT("P2P sending U from %d,%d down to %d,%d\n", thisIndex.x, thisIndex.y, i,thisIndex.y);
//	 thisProxy(i,thisIndex.y).updateRecvU(givenU);
//     }
    
  }
  
  //broadcast the L rightwards to the blocks in the same row
  inline void multicastRecvL() {
    traceUserSuppliedData(internalStep);
    traceMemoryUsage();
    
    DEBUG_PRINT("[PE %d] elem %d,%d Multicast to part of row %d from step %d\n", CkMyPe(), thisIndex.x, thisIndex.y, thisIndex.x, internalStep);
    
    CProxySection_LUBlk oneRow = CProxySection_LUBlk::ckNew(thisArrayID, thisIndex.x, thisIndex.x, 1, thisIndex.y+1, numBlks-1, 1);
    
    if (whichMulticastStrategy > -1)
      ComlibAssociateProxy(multicastStats[whichMulticastStrategy], oneRow);
    
    blkMsg *givenL = createABlkMsg();
    *(int*)CkPriorityPtr(givenL) = -1;
    oneRow.recvL(givenL);
    
//     for(int i=thisIndex.y+1; i<numBlks; i++){
//	 blkMsg *givenL = createABlkMsg();
//	 DEBUG_PRINT("P2P sending L from %d,%d right to %d,%d\n", thisIndex.x, thisIndex.y, thisIndex.x, i);
//	 thisProxy(thisIndex.x, i).updateRecvL(givenL);
//     }
    
  }

  void processComputeU(int ignoredParam) {
    DEBUG_PRINT("processComputeU() called on block %d,%d\n", thisIndex.x, thisIndex.y);
    CkAssert(internalStep==thisIndex.x && L);
    // We are in the top row of active blocks, and we
    // have received the incoming L
	
    // CkPrintf("[%d] chare %d,%d internalStep=%d computeU\n", CkMyPe(), thisIndex.x, thisIndex.y, internalStep);
    computeU(L);
    
    DEBUG_PRINT("[%d] chare %d,%d is top block this step, multicast U\n", CkMyPe(), thisIndex.x, thisIndex.y);
    multicastRecvU(); //broadcast the newly computed U downwards to the blocks in the same column
    
    CmiFree(UsrToEnv(L));
    
    DEBUG_PRINT("chare %d,%d is now done\n",  thisIndex.x, thisIndex.y);
  }

  void localSolve(double *xvec, double *preVec) {
      for (int i = 0; i < BLKSIZE; i++) {
        xvec[i] = 0.0;
      }
      
      for (int i = 0; i < BLKSIZE; i++) {
        for (int j = 0; j < BLKSIZE; j++) {
          xvec[i] += LU[i][j] * preVec[j];
        }
      }
  }

  void localForward(double *xvec) {
    for (int i = 0; i < BLKSIZE; i++) {
      for (int j = 0; j < i; j++) {
	xvec[i] -= LU[i][j] * xvec[j];
      }
    }
  }

  void localBackward(double *xvec) {
    for (int i = BLKSIZE-1; i >= 0; i--) {
      for (int j = i+1; j < BLKSIZE; j++) {
	xvec[i] -= LU[i][j] * xvec[j];
      }
      xvec[i] /= LU[i][i];
    }
  }


  void beginForward(int size, double *preVec) {
      // Perform local solve and reduce left-of-diagonal row to diagonal
      double *xvec = new double[BLKSIZE];
      localSolve(xvec, preVec);
      CkCallback cb(CkIndex_LUBlk::recvSolveData(0), thisProxy(thisIndex.x, thisIndex.x));
      mcastMgr->contribute(sizeof(double) * BLKSIZE, xvec, CkReduction::sum_double, rowBeforeCookie, cb, thisIndex.x);
  }


  void beginBackward(int size, double *preVec) {
      // Perform local solve and reduce right-of-diagonal row to diagonal
      double *xvec = new double[BLKSIZE];
      localSolve(xvec, preVec);
      CkCallback cb(CkIndex_LUBlk::recvSolveData(0), thisProxy(thisIndex.x, thisIndex.x));
      mcastMgr->contribute(sizeof(double) * BLKSIZE, xvec, CkReduction::sum_double, rowAfterCookie, cb, thisIndex.x);
  }


  void print() {
    this->print("LU-solution");
  }

  void print(const char* step) {
#if defined(PRINT_VECTORS)
    char buf[200];
    sprintf(buf, "%s-%d-%d", step, thisIndex.x, thisIndex.y);

    FILE *file = fopen(buf, "w+");

    for (int i = 0; i < BLKSIZE; i++) {
      for (int j = 0; j < BLKSIZE; j++) {
	fprintf(file, "%f ", LU[i][j]);
      }
      fprintf(file, "\n");
    }

    fclose(file);
#endif
  }

private:

  // Copy received pivot data into its place in this block
  void applySwap(int row, int offset, double *data, double b) {
      DEBUG_PIVOT("(%d, %d): remote pivot inserted at %d\n", thisIndex.x, thisIndex.y, row);
      bvec[row] = b;
      memcpy( &(LU[row][offset]), data, sizeof(double)*(BLKSIZE-offset) );
  }

  // Exchange local data
  void swapLocal(int row1, int row2, int offset=0) {
    if (row1 == row2) return;
    std::swap(bvec[row1], bvec[row2]);
    /// @todo: Is this better or is it better to do 3 memcpys
    std::swap_ranges( &(LU[row1][offset]), &(LU[row1][BLKSIZE]), &(LU[row2][offset]) );
  }

  void doPivotLocal(int row1, int row2) {
    // The chare indices where the two rows are located
    row1Index = row1 / BLKSIZE;
    row2Index = row2 / BLKSIZE;
    // The local indices of the two rows within their blocks
    localRow1 = row1 % BLKSIZE;
    localRow2 = row2 % BLKSIZE;
    remoteSwap = false;

    // If I hold portions of both the current row and pivot row, its a local swap
    if (row1Index == thisIndex.x && row2Index == thisIndex.x) {
      swapLocal(localRow1, localRow2);
      // else if I hold portions of at just one row, its a remote swap
    } else if (row1Index == thisIndex.x) {
      thisLocalRow = localRow1;
      otherRowIndex = row2Index;
      globalThisRow = row1;
      globalOtherRow = row2;
      remoteSwap = true;
    } else if (row2Index == thisIndex.x) {
      thisLocalRow = localRow2;
      otherRowIndex = row1Index;
      globalThisRow = row2;
      globalOtherRow = row1;
      remoteSwap = true;
    }
    // else, I dont have any data affected by this pivot op
  }


  /// Record the effect of a pivot operation in terms of actual row numbers
  void recordPivot(const int r1, const int r2)
  {
      numRowsSinceLastPivotSend++;
      // If the two rows are the same, then dont record the pivot operation at all
      if (r1 == r2) return;
      std::map<int,int>::iterator itr1, itr2;
      // The records for the two rows (already existing or freshly created)
      itr1 = (pivotRecords.insert( std::make_pair(r1,r1) ) ).first;
      itr2 = (pivotRecords.insert( std::make_pair(r2,r2) ) ).first;
      // Swap the values (the actual rows living in these two positions)
      std::swap(itr1->second, itr2->second);
  }



  /** Is it time to send out the next batch of pivots
   *
   * @note: Any runtime adaptivity should be plugged here
   */
  bool shouldSendPivots()
  {
      return (numRowsSinceLastPivotSend >= suggestedPivotBatchSize);
  }



  /// Periodically send out the agglomerated pivot operations
  void announceAgglomeratedPivots()
  {
      #if defined(VERBOSE_PIVOT_RECORDING) || defined(VERBOSE_PIVOT_AGGLOM)
          std::stringstream pivotLog;
          pivotLog<<"["<<thisIndex.x<<","<<thisIndex.y<<"]"
                  <<" announcing "<<pivotRecords.size()<<" pivot operations in batch "<<pivotBatchTag<<std::endl;
      #endif

      // Create and initialize a msg to carry the pivot sequences
      pivotSequencesMsg *msg = new(numRowsSinceLastPivotSend+1, numRowsSinceLastPivotSend*2, sizeof(int)*8)
                                  pivotSequencesMsg(pivotBatchTag, numRowsSinceLastPivotSend);
      msg->numSequences = 0;
      memset(msg->seqIndex,      0, sizeof(int) * numRowsSinceLastPivotSend );
      memset(msg->pivotSequence, 0, sizeof(int) * numRowsSinceLastPivotSend*2 );

      /// Parse the pivot operations and construct optimized pivot sequences
      int seqNo =-1, i = 0;
      std::map<int,int>::iterator itr = pivotRecords.begin();
      while (itr != pivotRecords.end())
      {
          #ifdef VERBOSE_PIVOT_RECORDING
            pivotLog<<std::endl;
          #endif
          msg->seqIndex[++seqNo] = i;
          int chainStart = itr->first;
          msg->pivotSequence[i++] = chainStart;
          #ifdef VERBOSE_PIVOT_RECORDING
            pivotLog<<chainStart;
          #endif
          while (itr->second != chainStart)
          {
              msg->pivotSequence[i] = itr->second;
              #ifdef VERBOSE_PIVOT_RECORDING
                  pivotLog<<" <-- "<<itr->second;
              #endif
              std::map<int,int>::iterator prev = itr;
              itr = pivotRecords.find(itr->second);
              pivotRecords.erase(prev);
              i++;
          }
          pivotRecords.erase(itr);
          itr = pivotRecords.begin();
      }
      msg->seqIndex[++seqNo] = i; ///< @note: Just so that we know where the last sequence ends
      msg->numSequences = seqNo;
      #if defined(VERBOSE_PIVOT_RECORDING) || defined(VERBOSE_PIVOT_AGGLOM)
          CkPrintf("%s\n", pivotLog.str().c_str());
      #endif

      // Make a copy of the msg for the left section too
      pivotSequencesMsg *leftMsg = (pivotSequencesMsg*) CkCopyMsg((void**)&msg);

      // Send the pivot ops to the right section (trailing sub-matrix chares + post-diagonal active row chares)
      *(int*)CkPriorityPtr(msg) = (thisIndex.y + 1) * BLKSIZE;
      CkSetQueueing(msg, CK_QUEUEING_IFIFO);
      pivotRightSection.applyPivots(msg);
      // Send the pivot ops to the left section (left of activeColumn and below activeRow)
      *(int*)CkPriorityPtr(leftMsg) = BLKSIZE * numBlks + 1;
      CkSetQueueing(leftMsg, CK_QUEUEING_IFIFO);
      pivotLeftSection.applyPivots(leftMsg);

      // Prepare for the next batch of agglomeration
      pivotRecords.clear();
      pivotBatchTag += numRowsSinceLastPivotSend;
      numRowsSinceLastPivotSend = 0;
  }


  /// Given a set of pivot ops, send out participating row chunks that you own
  void sendPendingPivots(pivotSequencesMsg *msg)
  {
      #ifdef VERBOSE_PIVOT_AGGLOM
          std::stringstream pivotLog;
          pivotLog<<"["<<thisIndex.x<<","<<thisIndex.y<<"]"
                  <<" processing "<<msg->numSequences
                  <<" pivot sequences in batch "<<pivotBatchTag;
      #endif

      int *pivotSequence = msg->pivotSequence, *idx = msg->seqIndex;
      int numSequences = msg->numSequences;

      // Count the number of rows that I send to each chare
      int numMsgsTo[numBlks];
      memset(numMsgsTo, 0, sizeof(int)*numBlks);
      for (int i=0; i< numSequences; i++)
      {
          for (int j=idx[i]; j<idx[i+1]; j++)
          {
              int recverIdx = pivotSequence[j]/BLKSIZE;
              int senderIdx =-1;
              // circular traversal of pivot sequence
              if (j < idx[i+1]-1)
                  senderIdx = pivotSequence[j+1]/BLKSIZE;
              else
                  senderIdx = pivotSequence[idx[i]]/BLKSIZE;
              // If I am the sending to another chare
              if (thisIndex.x == senderIdx && thisIndex.x != recverIdx)
                  numMsgsTo[recverIdx]++;
          }
      }

      // Preallocate msgs that are big enough to carry the rows I'll be sending to each of the other chares
      pivotRowsMsg* outgoingPivotMsgs[numBlks];
      memset(outgoingPivotMsgs, 0, sizeof(pivotRowsMsg*) * numBlks);
      for (int i=0; i< numBlks; i++)
      {
          // If this other chare expects row data from me
          if ( numMsgsTo[i] > 0)
          {
              // Create a big enough msg to carry all the rows I'll be sending to this chare
              outgoingPivotMsgs[i] = new(numMsgsTo[i], numMsgsTo[i]*BLKSIZE, numMsgsTo[i], sizeof(int)*8)
                                        pivotRowsMsg(BLKSIZE, pivotBatchTag);
              // Set a priority thats a function of your location wrt to the critical path
              *(int*)CkPriorityPtr(msg) = (thisIndex.y<internalStep)? numBlks*BLKSIZE : thisIndex.y * BLKSIZE;
              CkSetQueueing(msg, CK_QUEUEING_IFIFO);
          }
      }

      pendingIncomingPivots = 0;
      double *tmpBuf = NULL, tmpB;

      // Parse each sequence independently
      for (int i=0; i < numSequences; i++)
      {
          #ifdef VERBOSE_PIVOT_AGGLOM
              pivotLog<<"\n["<<thisIndex.x<<","<<thisIndex.y<<"] sequence "<<i<<": ";
          #endif

          // Find the location of this sequence in the msg buffer
          int *first      = pivotSequence + idx[i];
          int *beyondLast = pivotSequence + idx[i+1];
          CkAssert(beyondLast - first >= 2);

          // Identify a remote row in the pivot sequence as a point at which to
          // start and stop processing the circular pivot sequence
          int *ringStart = first;
          while ( (*ringStart / BLKSIZE == thisIndex.x) && (ringStart < beyondLast) )
              ringStart++;
          int *ringStop = ringStart;

          // If there are no remote rows in the sequence, we *have* to use a tmp buffer
          // The tmp buffer will now complete the circular sequence
          bool isSequenceLocal = false;
          if (ringStart == beyondLast)
          {
              isSequenceLocal = true;
              ringStart = first;
              ringStop  = beyondLast - 1;

              if (NULL == tmpBuf) tmpBuf = new double[BLKSIZE];
              memcpy(tmpBuf, LU[*ringStart%BLKSIZE], BLKSIZE*sizeof(double));
              tmpB = bvec[*ringStart%BLKSIZE];
              #ifdef VERBOSE_PIVOT_AGGLOM
                  pivotLog<<"tmp <-cpy- "<<*ringStart<<"; ";
              #endif
          }

          // Process all the pivot operations in the circular sequence
          int *to = ringStart;
          do
          {
              int *from        = (to+1 == beyondLast) ? first : to + 1;
              int fromChareIdx = *from / BLKSIZE;
              // If the current source row in the pivot sequence belongs to me, send it
              if (fromChareIdx == thisIndex.x)
              {
                  int toChareIdx = *to / BLKSIZE;
                  int fromLocal  = *from % BLKSIZE;
                  // If you're sending to yourself, memcopy
                  if (toChareIdx == thisIndex.x)
                  {
                      applySwap(*to%BLKSIZE, 0, LU[fromLocal], bvec[fromLocal]);
                      #ifdef VERBOSE_PIVOT_AGGLOM
                          pivotLog<<*to<<" <-cpy- "<<*from<<"; ";
                      #endif
                  }
                  // else, copy the data into the appropriate msg
                  else
                  {
                      outgoingPivotMsgs[*to/BLKSIZE]->copyRow(*to, LU[fromLocal], bvec[fromLocal]);
                      #ifdef VERBOSE_PIVOT_AGGLOM
                          pivotLog<<*to<<" <-msg- "<<*from<<"; ";
                      #endif
                  }
              }
              // else, the source data is remote
              else
              {
                  // if the current destination row belongs to me, make sure I expect the remote data
                  if (*to / BLKSIZE == thisIndex.x)
                  {
                      pendingIncomingPivots++;
                      #ifdef VERBOSE_PIVOT_AGGLOM
                          pivotLog<<*to<<" <-inwd- "<<*from<<"; ";
                      #endif
                  }
                  // else, i dont worry about this portion of the exchange sequence which is completely remote
                  else
                  {
                      #ifdef VERBOSE_PIVOT_AGGLOM
                          pivotLog<<*to<<" <-noop- "<<*from<<"; ";
                      #endif
                  }
              }
              // Setup a circular traversal of the pivot sequence
              if (++to == beyondLast)
                  to = first;
          } while (to != ringStop); // Keep going till you complete the ring

          // If the sequence was completely local, complete the circular sequence
          // by copying the temp buffer back into the matrix block
          if (isSequenceLocal)
          {
              applySwap(*(beyondLast-1)%BLKSIZE, 0, tmpBuf, tmpB);
              #ifdef VERBOSE_PIVOT_AGGLOM
                  pivotLog<<*(beyondLast-1)<<"<-cpy- tmp "<<"; ";
              #endif
          }

      } // end for loop through all sequences

      // Send out all the msgs carrying pivot data to other chares
      for (int i=0; i< numBlks; i++)
          if (numMsgsTo[i] > 0)
              thisProxy(i, thisIndex.y).acceptPivotData(outgoingPivotMsgs[i]);

      if (tmpBuf) delete [] tmpBuf;
      #ifdef VERBOSE_PIVOT_AGGLOM
          CkPrintf("%s\n",pivotLog.str().c_str());
      #endif
  }


  //internal functions for creating messages to encapsulate the priority
  inline blkMsg* createABlkMsg() {
    blkMsg *msg = mgr->createBlockMessage(thisIndex.x, thisIndex.y,
                                          internalStep, sizeof(int)*8);
    msg->setMsgData(LU[0], internalStep, BLKSIZE);
    return msg;
  }



  locval findLocVal(int startRow, int col, locval first = locval()) {
    locval l = first;
    for (int row = startRow; row < BLKSIZE; row++)
      if ( fabs(LU[row][col]) > fabs(l.val) ) {
        l.val = LU[row][col];
        l.loc = row + BLKSIZE * thisIndex.x;
      }
    return l;
  }



  /// Compute the multipliers and update the trailing sub-block
  void diagonalUpdate(int col) {
      // Only diagonal chares should use this method
      CkAssert(thisIndex.x == thisIndex.y);
    // Pivoting is done, so the diagonal entry better not be zero; else the matrix is singular
    if (fabs(LU[col][col]) <= 100 * std::numeric_limits<double>::epsilon() )
        CkAbort("Diagonal element very small despite pivoting. Is the matrix singular??");

    for(int j=col+1; j<BLKSIZE; j++) {
        LU[j][col] = LU[j][col]/ LU[col][col];
        for(int k=col+1;k<BLKSIZE;k++)
            LU[j][k] -= LU[j][col] * LU[col][k];
    }
  }



  /// Update the sub-block of this L block starting at specified
  /// offset from the active column
  void updateLsubBlock(int activeCol, double* U, int offset=1) {
      // Should only get called on L blocks
      CkAssert(thisIndex.x > thisIndex.y);
#if 0
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		  BLKSIZE, BLKSIZE-(activeCol+offset), 1,
		  -1.0, &LU[0][activeCol], BLKSIZE,
		  U+offset, BLKSIZE,
		  1.0, &LU[0][activeCol+offset], BLKSIZE);
#else
      for(int j = 0; j < BLKSIZE; j++)
          for(int k = activeCol+offset; k<BLKSIZE; k++)
              LU[j][k] -=  LU[j][activeCol] * U[k-activeCol];
#endif
  }


  /// Compute the multipliers based on the pivot value in the
  /// received row of U and also find the candidate pivot in
  /// the immediate next column (after updating it simultaneously)
  locval computeMultipliersAndFindColMax(int col, double *U)
  {
      // Should only get called on L blocks
      CkAssert(thisIndex.x > thisIndex.y);
      locval maxVal;

      if (col < BLKSIZE -1) {
          for (int j = 0; j < BLKSIZE; j++) {
              // Compute the multiplier
              LU[j][col]    = LU[j][col] / U[0];
              // Update the immediate next column
              LU[j][col+1] -= LU[j][col] * U[1];
              // Update the max value thus far
              if ( fabs(LU[j][col+1]) > fabs(maxVal.val) ) {
                  maxVal.val = LU[j][col+1];
                  maxVal.loc = j;
              }
          }
          // Convert local row num to global rownum
          maxVal.loc += thisIndex.x*BLKSIZE;
      }
      else {
          for (int j = 0; j < BLKSIZE; j++)
              LU[j][col]    = LU[j][col] / U[0];
      }

      return maxVal;
  }

};

#include "lu.def.h"
