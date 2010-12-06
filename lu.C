/**
   @file A Charm++ implementation of LU 
   
   Messages of type blkMsg contain blocks of data that are sent down and to the 
   right as the computation progresses. 

   As the blkMsg messages arrive, they are always put into a buffer, and a call 
   to thisProxy(thisIndex.x,thisIndex.y).progress() is made if the computation 
   can proceed without any additional messages.

 */

//#include <assert.h>
//#include <stdlib.h>
#include <string.h>
#include <iostream>
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

#elif USE_GSL_H
#include <gsl/gsl_cblas.h>

#else
#error "No BLAS Header files included!"
#endif

#if USE_ESSL_TEST
#define _ESVCPTR
#include <complex>
#include <essl.h>
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
#else
    #define DEBUG_PRINT(...)
    #define DEBUG_PIVOT(...)
#endif

#include <cmath>

class locval {
    public:
        locval(): val(0.0), loc(-1) {}
        locval(double _val, int _loc): val(_val), loc(_loc) {}
        double val;
        int loc;
};

void doExit(void *) { CkExit(); }

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

struct UMsg : public CMessage_UMsg, public CkMcastBaseMsg {
  double *data;
};

struct pivotMsg: public CMessage_pivotMsg, public CkMcastBaseMsg {
  int activeRow, row1, row2;

  pivotMsg(int activeRow_, int row1_, int row2_) :
    activeRow(activeRow_), row1(row1_), row2(row2_) {
    CkSetRefNum(this, activeRow);
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

static inline void takeRef(blkMsg *m) {
    CmiLock(lock);
    CmiReference(UsrToEnv(m));
    CmiUnlock(lock);
}
static inline void dropRef(blkMsg *m) {
    CmiLock(lock);
    CmiFree(UsrToEnv(m));
    CmiUnlock(lock);
}

class Main : public CBase_Main {
  double startTime;
  int iteration;
  int numIterations;
  int numBlks;

  int gMatSize;
  int BLKSIZE;
  int whichMulticastStrategy;
  int mapping;
  int memThreshold;
  bool solved, LUcomplete, workStarted;
  bool sentVectorData;

  CProxy_LUBlk luArrProxy;

public:
    Main(CkArgMsg* m) : iteration(0), numIterations(1), solved(false), LUcomplete(false), workStarted(false), sentVectorData(false) {

    if (m->argc<4) {
      CkPrintf("Usage: %s <matrix size> <block size> <mem threshold> [<iterations>]\n", m->argv[0]);
      CkExit();
    }

    gMatSize = atoi(m->argv[1]);
    BLKSIZE = atoi(m->argv[2]);
    memThreshold = atoi(m->argv[3]);

    if (m->argc > 3) {
      /*sscanf( m->argv[2], "%d", &strategy);
	CkPrintf("CLI: strategy=%d\n", strategy);*/
      
      /*
      sscanf( m->argv[2], "%d", &memThreshold);
      CkPrintf("CLI: memThreshold=%dMB\n", memThreshold);
      */


      if (m->argc >= 5)
	numIterations = atoi(m->argv[4]);
      CkPrintf("CLI: numIterations=%d\n", numIterations);
    }

    if (gMatSize%BLKSIZE!=0) {
      CkPrintf("The matrix size %d should be a multiple of block size %d!\n", gMatSize, BLKSIZE);
      CkExit();
    }
    numBlks = gMatSize / BLKSIZE;

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

    CkPrintf("Running LU on %d processors (%d nodes) on matrix %dX%d with control points\n",
	     CkNumPes(), CmiNumNodes(), gMatSize, gMatSize);

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
    CkStartQD(CkCallback(doExit));
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
    //CkCallback *cb = new CkCallback(CkIndex_Main::arrayIsCreated(NULL), thisProxy);
    //CkStartQD(*cb); // required currently for use with new Comlib
    thisProxy.arrayIsCreated();
  }

  void arrayIsCreated() {
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
    
    
      mapping = 2;

      char note[200];
      sprintf(note, "*** New iteration: block size = %d, mapping = %d %s, multicast = %d, memthreshold = %d MB", 
	      BLKSIZE, mapping, mapping == 1 ? "Balanced Snake" : "Block Cylic", whichMulticastStrategy, memThreshold);
      traceUserSuppliedNote(note);
      CkPrintf("%s\n", note);
      fflush(stdout);
    
      CkArrayOptions opts(numBlks, numBlks);
      opts.setAnytimeMigration(false)
	  .setStaticInsertion(true);
      switch (mapping) {
      case 0:
	opts.setMap(CProxy_BlockCyclicMap::ckNew());
	break;
      case 1:
        opts.setMap(CProxy_LUBalancedSnakeMap2::ckNew(numBlks, BLKSIZE));
	break;
      case 2:
	  opts.setMap(CProxy_RealBlockCyclicMap::ckNew(1, numBlks));
	  break;
      }
      CProxy_LUMgr mgr = CProxy_PrioLU::ckNew(BLKSIZE, gMatSize);

      luArrProxy = CProxy_LUBlk::ckNew(opts);

      workStarted = true;
    
      luArrProxy.startup(0, BLKSIZE, numBlks, memThreshold, mgr);
    }
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
    std::cout << "RESULT procs: \t" << CkNumPes() << "\tblock size:\t" << BLKSIZE << "\tGFlops:\t" << gflops << "\tTime(s):\t" << duration << std::endl;

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

	double n = BLKSIZE * numBlks;

	double r = maxvals[3]/((maxvals[0]*maxvals[2]+maxvals[1])*n*std::numeric_limits<double>::epsilon());

	DEBUG_PRINT("|A|inf = %e\n|b|inf = %e\n|x|inf = %e\n|Ax-b|inf = %e\n",maxvals[0],maxvals[1],maxvals[2],maxvals[3]);

	CkPrintf("epsilon = %e\nr = %f\n",std::numeric_limits<double>::epsilon(),r);
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
  // The section of chares in the array on and below the current diagonal
  CProxySection_LUBlk belowLeft, belowRight;

  /// Variables used during factorization
  double *LU;

  int BLKSIZE, numBlks;
  blkMsg *L, *U;
  int internalStep, activeCol, currentStep, ind;

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
  locval l;
  int pivotBlk;

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

	  //Regenerate A and place into already allocated LU
    genBlock();

	  //Diagonals regenerate b and distribute x across entire column
	  if(thisIndex.x == thisIndex.y)
	  {
		  b = new double[BLKSIZE];
		  genVec(b);

		  CProxySection_LUBlk col =
				  CProxySection_LUBlk::ckNew(thisArrayID, 0, numBlks-1,
						  1, thisIndex.y, thisIndex.y, 1);

          col.recvXvec(BLKSIZE, bvec);
          Ax = new double[BLKSIZE];

	  }
  }

  //VALIDATION
  void recvXvec(int size, double* xvec) {

	  double *partial_b = new double[BLKSIZE];

	  //Perform local dgemv
#if USE_ESSL
    dgemv("T", BLKSIZE, BLKSIZE, 1.0, LU, BLKSIZE, xvec, 1, 0.0, partial_b, 1);
#else
    cblas_dgemv( CblasRowMajor, CblasNoTrans,
    		  BLKSIZE, BLKSIZE, 1.0, LU,
    		  BLKSIZE, xvec, 1, 0.0, partial_b, 1);
#endif

    //sum-reduction of result across row with diagonal element as target
    thisProxy(thisIndex.x,thisIndex.x).sumBvec(BLKSIZE,partial_b);

    //if you are not the diagonal, find your max A value and contribute
    if(thisIndex.x != thisIndex.y) {
    	//find local max of A
    	double A_max = infNorm(BLKSIZE * BLKSIZE, LU);
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
		  for(int i = 0; i < BLKSIZE; i++)
			  Ax[i] = 0.0;
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
	  double *residuals = new double[BLKSIZE];

	  //diagonal elements that received sum-reduction perform b - A*x
	  for (int i = 0; i < BLKSIZE; i++) {
		  residuals[i] = b[i] - Ax[i];
//		  if(fabs(residuals[i]) > 1e-14 || std::isnan(residuals[i]) || std::isinf(residuals[i]))
//			  CkPrintf("WARNING: Large Residual for x[%d]: %f - %f = %e\n", thisIndex.x*BLKSIZE+i, b[i], bvec[i], residuals[i]);
	  }

	  //find local max values
	  double A_max = infNorm(BLKSIZE * BLKSIZE, LU);
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

  void genBlock()
    {
      MatGen rnd(seed_A);

      // Skip the rows of blocks before this one
      rnd.skipNDoubles(thisIndex.x * BLKSIZE * BLKSIZE * numBlks);
      for (int row = 0; row < BLKSIZE; ++row) {
	// Skip the blocks before this one in this row
	rnd.skipNDoubles(thisIndex.y * BLKSIZE);
	rnd.getNRndDoubles(BLKSIZE, LU + row*BLKSIZE);
	// Skip the blocks after this one in this row
	rnd.skipNDoubles((numBlks - thisIndex.y - 1) * BLKSIZE);
      }

    }

  void genVec(double *buf)
    {
      MatGen rnd(seed_b);

      // Skip the blocks before this one
      rnd.skipNDoubles(thisIndex.x * BLKSIZE);
      rnd.getNRndDoubles(BLKSIZE, buf);
    }

    void init(int _whichMulticastStrategy, int _BLKSIZE, int _numBlks,
	      int memThreshold, CProxy_LUMgr _mgr) {
    whichMulticastStrategy = _whichMulticastStrategy;
    BLKSIZE = _BLKSIZE;
    numBlks = _numBlks;
    mgr = _mgr.ckLocalBranch();
    
    // Set the schedulers memory usage threshold to the one based upon a control point
    schedAdaptMemThresholdMB = memThreshold;

    CkAssert(BLKSIZE>0); // If this fails, readonly variables aren't
			 // propagated soon enough. I'm assuming they
			 // are safe to use here.

#if USE_MEMALIGN
    LU = (double*)memalign(128, BLKSIZE*BLKSIZE*sizeof(double) );
    //	 CkPrintf("LU mod 128 = %lu\n", ((unsigned long)LU) % 128);
    CkAssert(LU != NULL);
#else
    LU = new double[BLKSIZE*BLKSIZE];
#endif

    internalStep = 0;  
     
    traceUserSuppliedData(-1);	
    traceMemoryUsage();	 
     
    //VALIDATION: saved seed value to use for validation
    seed_A = 2998389;
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
        activePanel.prepareForPivotRedn(activePanelMsg);
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
    free(LU);
#else
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

#if USE_ESSL_TEST
    dtrsm("R", "U", "N", "U", BLKSIZE, BLKSIZE, 1.0, givenL, BLKSIZE, LU, BLKSIZE);
#else
    cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, BLKSIZE, BLKSIZE, 1.0, givenL, BLKSIZE, LU, BLKSIZE);
#endif
  }

  void computeL(blkMsg *givenUMsg) {
    traceLU t(internalStep, traceComputeL);
    double *givenU = givenUMsg->data;
    //	  CkAssert( ((unsigned long)givenU) % 16 == 0);

    DEBUG_PRINT("elem[%d,%d]::computeL called at step %d\n", thisIndex.x, thisIndex.y, internalStep);

#if USE_ESSL
    dtrsm("R", "U", "N", "N", BLKSIZE, BLKSIZE, 1.0, givenU, BLKSIZE, LU, BLKSIZE);
#else
    cblas_dtrsm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, BLKSIZE, BLKSIZE, 1.0, givenU, BLKSIZE, LU, BLKSIZE);
#endif
  }

  void updateMatrix(blkMsg *givenLMsg, blkMsg *givenUMsg) {
    traceLU t(internalStep, traceTrailingUpdate);
    DEBUG_PRINT("elem[%d,%d] is updating its value at step %d\n", thisIndex.x, thisIndex.y, internalStep);

    double *incomingL = givenLMsg->data;
    double *incomingU = givenUMsg->data;

#if USE_ESSL
    dgemm( "N", "N",
	   BLKSIZE, BLKSIZE, BLKSIZE,
	   -1.0, incomingU,
	   BLKSIZE, incomingL, BLKSIZE,
	   1.0, LU, BLKSIZE);
#else
    cblas_dgemm( CblasRowMajor,
		 CblasNoTrans, CblasNoTrans,
		 BLKSIZE, BLKSIZE, BLKSIZE,
		 -1.0, incomingL,
		 BLKSIZE, incomingU, BLKSIZE,
		 1.0, LU, BLKSIZE);
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

  int inline getIndex(int i, int j) {
    return i * BLKSIZE + j;
  }

  void localSolve(double *xvec, double *preVec) {
      for (int i = 0; i < BLKSIZE; i++) {
        xvec[i] = 0.0;
      }
      
      for (int i = 0; i < BLKSIZE; i++) {
        for (int j = 0; j < BLKSIZE; j++) {
          xvec[i] += LU[getIndex(i,j)] * preVec[j];
        }
      }
  }

  void localForward(double *xvec) {
    for (int i = 0; i < BLKSIZE; i++) {
      for (int j = 0; j < i; j++) {
	xvec[i] -= LU[getIndex(i,j)] * xvec[j];
      }
    }
  }

  void localBackward(double *xvec) {
    for (int i = BLKSIZE-1; i >= 0; i--) {
      for (int j = i+1; j < BLKSIZE; j++) {
	xvec[i] -= LU[getIndex(i,j)] * xvec[j];
      }
      xvec[i] /= LU[getIndex(i,i)];
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
	fprintf(file, "%f ", LU[getIndex(i,j)]);
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
    for (int col = offset; col < BLKSIZE; ++col)
      LU[getIndex(row, col)] = data[col - offset];
  }

  // Exchange local data
  void swapLocal(int row1, int row2, int offset=0) {
    double buf;
    buf = bvec[row1];
    bvec[row1] = bvec[row2];
    bvec[row2] = buf;
    // Swap the row of A (LU)
    for (int col = offset; col < BLKSIZE; col++) {
      buf = LU[getIndex(row1, col)];
      LU[getIndex(row1, col)] = LU[getIndex(row2, col)];
      LU[getIndex(row2, col)] = buf;
    }
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

  //internal functions for creating messages to encapsulate the priority
  inline blkMsg* createABlkMsg() {
    blkMsg *msg = mgr->createBlockMessage(thisIndex.x, thisIndex.y,
                                          internalStep, sizeof(int)*8);
    msg->setMsgData(LU, internalStep, BLKSIZE);
    return msg;
  }

  locval findLocVal(int startRow, int col, locval first = locval()) {
    locval l = first;
    for (int row = startRow; row < BLKSIZE; row++)
      if ( fabs(LU[getIndex(row, col)]) > fabs(l.val) ) {
        l.val = LU[getIndex(row, col)];
        l.loc = row + BLKSIZE * thisIndex.x;
      }
    return l;
  }

  // Local multiplier computation and update after U is sent to the blocks below
  void diagonalUpdate(int col) {
    // Pivoting is done, so the diagonal entry better not be zero; else the matrix is singular
    if (fabs(LU[getIndex(col,col)]) <= 100 * std::numeric_limits<double>::epsilon() )
        CkAbort("Diagonal element very small despite pivoting. Is the matrix singular??");

    computeMultipliers(LU[getIndex(col,col)],col+1,col);
    for(int k=col+1;k<BLKSIZE;k++) {
      for(int j=col+1; j<BLKSIZE; j++) {
        LU[getIndex(j,k)] = LU[getIndex(j,k)] - LU[getIndex(j,col)] * LU[getIndex(col, k)];
      }
    }
  }

  /// Compute the multipliers based on the pivot value from the diagonal chare
  //  starting at [row, col]
  void computeMultipliers(double a_kk, int row, int col) {
#if 0
      cblas_dscal(BLKSIZE-row, 1/a_kk, &LU[getIndex(row,col)], BLKSIZE);
#else
    for(int i = row; i<BLKSIZE;i++)
      LU[getIndex(i,col)] = LU[getIndex(i,col)]/a_kk;
#endif
  }

  /// Update the values in the columns ahead of the active column within the
  // pivot section based on the Usegment and the multipliers
  void updateAllCols(int col, double* U) {
    //TODO: Replace with DGEMM
#if 0
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		  BLKSIZE, BLKSIZE-(col+1), 1,
		  -1.0, &LU[getIndex(0,col)], BLKSIZE,
		  U+1, BLKSIZE,
		  1.0, &LU[getIndex(0,col+1)], BLKSIZE);
#else
    for(int k=col+1;k<BLKSIZE;k++) {
      for(int j=0; j<BLKSIZE; j++) {
        //U[k] might need to be U[j]?
        LU[getIndex(j,k)] =  LU[getIndex(j,k)] - LU[getIndex(j,col)]*U[k-col];
      }
    }
#endif
  }

};

#include "lu.def.h"
