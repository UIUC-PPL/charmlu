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


#if USE_MEMALIGN
#include <malloc.h>
#endif

#include <comlib.h>
#include <controlPoints.h> // must come before user decl.h if they are using the pathInformationMsg
#include "lu.decl.h"
#include <trace-projections.h>

#include <queueing.h> // for access to memory threshold setting

/* readonly: */
CProxy_Main mainProxy;
CProxy_LUBlk luArrProxy;
int gMatSize;
int traceTrailingUpdate;
int traceComputeU;
int traceComputeL;
int traceSolveLocalLU;
int doPrioritize;
ComlibInstanceHandle multicastStats[4];
 
//#define DEBUG_PRINT(...) CkPrintf(__VA_ARGS__)
#define DEBUG_PRINT(...) 

struct locval {
  double val;
  int loc;
};

CkReductionMsg *maxLocVal(int nMsg, CkReductionMsg **msgs)
{
  CkAssert(nMsg > 0);

  double val;
  int loc;
  locval l = *((locval *)(msgs[0]->getData()));

  for (int i = 1; i < nMsg; ++i) {
    loc = ((locval *)(msgs[i]->getData()))->loc;
    val = ((locval *)(msgs[i]->getData()))->val;

    if (val > l.val) {
      l.val = val;
      l.loc = loc;
    }
  }

  return CkReductionMsg::buildNew(sizeof(locval), &l);
}

enum continueWithTask {
  NO_CONTINUE = 0,
  CONTINUE_LU,
  CONTINUE_U,
  CONTINUE_L,
  CONTINUE_TRAIL
};

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

  void getNRndInts(int num, int *d) {
    for (int i=0; i<num; i++)
      d[i] = nextRndInt();
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

/** do a space filling curve style allocation from the bottom right to the upper left. */
class LUSnakeMap: public CkArrayMap {
public:
  int numBlks, BLKSIZE;
  LUSnakeMap(int _numBlks, int _BLKSIZE) : numBlks(_numBlks), BLKSIZE(_BLKSIZE) {}
  int procNum(int arrayHdl, const CkArrayIndex &idx) {
    int *coor = (int *)idx.data();
    int x = coor[0];
    int y = coor[1];

    int numProcs = CkNumPes();

    const int numsteps=numBlks;
   
    int p=0;
    for(int step=numsteps-1;step>=0;step--){
      
      // go along row
      for(int i=1;i<numsteps-step;i++){
	int curr_x = i + step;
	int curr_y = step;
	if(x==curr_x && y==curr_y)
	  return p % numProcs;
	p++;
      }
      
      // visit corner
      if(x==step && y==step)
	return p % numProcs;
      p++;
      
      // go along column
      for(int i=1;i<numsteps-step;i++){
	int curr_y = i + step;
	int curr_x = step;
	if(x==curr_x && y==curr_y)
	  return p % numProcs;
	p++;
      }
    }
    
    CkAbort("Mapping code has a bug in it\n");
    return -1;
  }
};

/** do an allocation that results in almost identical numbers of trailing updates per PE */
class LUBalancedSnakeMap: public CkArrayMap {
public:
  
  int mappingSize;
  int *mapping;
  int *peLoads;
  int numBlks, BLKSIZE;
  

  void setMapping(int x, int y, int pe){
    CkAssert(y*numBlks+x < mappingSize);
    mapping[y*numBlks+x] = pe;
    peLoads[pe] += workLoad(x, y);
  }

  int getMapping(int x, int y){
    CkAssert(y*numBlks+x < mappingSize);
    return mapping[y*numBlks+x];
  }

  /** build and store the mapping once */
  LUBalancedSnakeMap(int _numBlks, int _BLKSIZE) : numBlks(_numBlks), BLKSIZE(_BLKSIZE) {
    int numProcs = CkNumPes();

    mappingSize = numBlks*numBlks;
    mapping = new int[mappingSize];    

    peLoads = new int[numProcs];
    for(int i=0;i<numProcs;i++)
      peLoads[i]=0;
    
    const int numsteps=numBlks;
    

    // go along first row
    int step = 0;
    for(int i=1;i<numsteps-step;i++){
      int x = i + step;
      int y = step;
      setMapping(x, y, minLoadedPE() );
    }
    
    // visit first corner & column
    for(int i=0;i<numsteps-step;i++){
      int y = i + step;
      int x = step;
      setMapping(x, y, minLoadedPE() );
    }
    
    for(int step=numsteps-1;step>=1;step--){
      
      // go along row
      for(int i=1;i<numsteps-step;i++){
	int x = i + step;
	int y = step;
	setMapping(x, y, minLoadedPE() );
      }
      
      // visit corner & column
      for(int i=0;i<numsteps-step;i++){
	int y = i + step;
	int x = step;
	setMapping(x, y, minLoadedPE() );
      }
      
    }
    
  }
  
  int minLoadedPE(){
    int minLoadFound = 1000000000;
    int minPEFound = CkNumPes()-1;
    for(int p=CkNumPes()-1; p>=0; p--){
      if(peLoads[p] < minLoadFound) {
	minPEFound = p;
	minLoadFound = peLoads[p];
      }
    }
    return minPEFound;	  
  }
  
  int workLoad(int x, int y){
    if(x<y)
      return x+1;
    else 
      return y+1;
  }
  
  int procNum(int arrayHdl, const CkArrayIndex &idx) {
    int *coor = (int *)idx.data();
    int x = coor[0];
    int y = coor[1];
    return getMapping(x,y);
  }

};

/** do an allocation that results in almost identical numbers of trailing updates per PE */
class LUBalancedSnakeMap2: public CkArrayMap {
public:
  int mappingSize;
  int *mapping;
  int *peLoads;
  int stateN;
  int numBlks, BLKSIZE;

  void setMapping(int x, int y, int pe){
    CkAssert(y*numBlks+x < mappingSize);
    mapping[y*numBlks+x] = pe;
    peLoads[pe] += workLoad(x, y);
  }

  int getMapping(int x, int y){
    CkAssert(y*numBlks+x < mappingSize);
    return mapping[y*numBlks+x];
  }

  /** build and store the mapping once */
  LUBalancedSnakeMap2(int _numBlks, int _BLKSIZE) : numBlks(_numBlks), BLKSIZE(_BLKSIZE) {
    stateN = 0;
    int numProcs = CkNumPes();

    mappingSize = numBlks*numBlks;
    mapping = new int[mappingSize];  

    peLoads = new int[numProcs];
    for (int i = 0; i < numProcs; i++)
      peLoads[i]=0;
    
    const int numsteps = numBlks;

    for(int step = numsteps - 1; step >= 2; step--){
      // go along row
      for(int i = 1; i < numsteps - step; i++){
	int x = i + step;
	int y = step;
	int minLoaded = minLoadedPE();
	setMapping(x, y, minLoaded);
      }
      
      // visit corner & column
      for(int i = 0; i < numsteps - step; i++){
	int y = i + step;
	int x = step;
	int minLoaded = minLoadedPE();
	setMapping(x, y, minLoaded);
      }
    }

    // go along first two rows
    for (int x = 1; x < numsteps; x++){
      int minLoaded = minLoadedPE();
      setMapping(x, 0, minLoaded);
      setMapping(x, 1, minLoaded);
    }
    
    // visit first corner & first two columns
    for (int y = 0; y < numsteps; y++) {
      int minLoaded = minLoadedPE();
      setMapping(0, y, minLoaded);
      if (y != 0)
	setMapping(1, y, minLoaded);
    }

  }
  
  int minLoadedPE() {
    int minLoadFound = 1000000000;
    int minPEFound = CkNumPes()-1;
    for(int p = CkNumPes() - 1; p >= 0; p--){
      if(peLoads[p] < minLoadFound) {
	minPEFound = p;
	minLoadFound = peLoads[p];
      }
    }

    if (stateN == minPEFound) {
      int proc = stateN++ % CkNumPes();
      stateN = proc;
      return proc;
    }

    stateN = minPEFound;
    return minPEFound;	  
  }
  
  int workLoad(int x, int y){
    if (x < y)
      return x+1;
    else 
      return y+1;
  }
  
  int procNum(int arrayHdl, const CkArrayIndex &idx) {
    int *coor = (int *)idx.data();
    int x = coor[0];
    int y = coor[1];
    return getMapping(x,y);
  }

};

class BlockCyclicMap: public CkArrayMap {
public:
  BlockCyclicMap() {}
  int procNum(int arrayHdl, const CkArrayIndex &idx) {
    int *coor = (int *)idx.data();
    int x = coor[0];
    int y = coor[1];

    int numProcs = CkNumPes();
    int numNodes = CkNumNodes();
    //int numNodes = numProcs/4;

    int procPerNode = numProcs/numNodes;
    //int procPerNode = 4;
    CkAssert(numProcs%numNodes==0);

    //assume a procPerNode X numNodes grid of processors

    int penum = (x%procPerNode)*numNodes+(y%numNodes);
    return penum;
  }
};

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

class Main : public CBase_Main {
public:
  double startTime;
  int iteration;
  int numIterations;
  int numBlks;

  int BLKSIZE;
  int whichMulticastStrategy;
  int mapping;
  int memThreshold;
  bool solved, LUcomplete, workStarted;
  bool sentVectorData;

  Main(CkArgMsg* m) : numIterations(1), solved(false), LUcomplete(false), workStarted(false), sentVectorData(false) {
    iteration = 0;

    if (m->argc<3) {
      CkPrintf("Usage: %s <matrix size> <mem threshold> <iterations>\n", m->argv[0]);
      CkExit();
    }

    if (m->argc > 3) {
      /*sscanf( m->argv[2], "%d", &strategy);
	CkPrintf("CLI: strategy=%d\n", strategy);*/
      
      /*
      sscanf( m->argv[2], "%d", &memThreshold);
      CkPrintf("CLI: memThreshold=%dMB\n", memThreshold);
      */


      if (m->argc >= 4)
	numIterations = atoi(m->argv[3]);
      CkPrintf("CLI: numIterations=%d\n", numIterations);
    }

    mainProxy = thisProxy;
      
    DEBUG_PRINT("Registering user events\n");
    traceTrailingUpdate = traceRegisterUserEvent("Trailing Update");
    traceComputeU = traceRegisterUserEvent("Compute U");
    traceComputeL = traceRegisterUserEvent("Compute L");
    traceSolveLocalLU = traceRegisterUserEvent("Solve local LU");
    
    traceRegisterUserEvent("Local Multicast Deliveries", 10000);    
    traceRegisterUserEvent("Remote Multicast Forwarding - preparing", 10001);
    traceRegisterUserEvent("Remote Multicast Forwarding - sends", 10002);
    
    gMatSize = atoi(m->argv[1]);
 
    /*if (gMatSize%1024!=0) 
      CkAbort("The matrix size should be a multiple of 1024!\n");*/
  
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


    thisProxy.iterationCompleted();

  }

  void finishInit() {
    if (!sentVectorData) {
      // Generate vector data
      double* vec = new double[BLKSIZE*numBlks];
    
      for (int i = 0; i < BLKSIZE*numBlks; i++)
	vec[i] = i;
      
      for (int i = 0; i < numBlks; i++)
        for (int j = 0; j < numBlks; j++)
          luArrProxy(i, j).initVec(BLKSIZE, vec+(i*BLKSIZE));

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

      //luArrProxy(0, 0).print();

      CkPrintf("called solve()\n");

      // Perform forward solving
      luArrProxy(0, 0).solve(false, 0, NULL);

      solved = true;
      iteration++;
    } else {
      // Only advance phases after a few factorizations have been performed
      // Prior to the first phase of actual work, iteration=1
      if( 1 || iteration % 2 == 1 || iteration==1){
	gotoNextPhase();
      
	whichMulticastStrategy = controlPoint("multicast_strategy", 2, 2);
	BLKSIZE = 1 << 2; //1 << controlPoint("block_size", 10,10);
	mapping = controlPoint("mapping", 1, 1);
	memThreshold = 200 + controlPoint("memory_threshold", 0, 20) * 100;
      
	// CkPrintf("%d %d %d\n",  (int)BLKSIZE, (int)mapping, (int)whichMulticastStrategy);
	// fflush(stdout);
	// fflush(stderr);
      
	numBlks = gMatSize/BLKSIZE;
      }
    
    
      char note[200];
      sprintf(note, "*** New iteration: block size = %d, mapping = %d %s, multicast = %d, memthreshold = %d MB", 
	      BLKSIZE, mapping, mapping == 1 ? "Balanced Snake" : "Block Cylic", whichMulticastStrategy, memThreshold);
      traceUserSuppliedNote(note);
      CkPrintf("%s\n", note);
      fflush(stdout);
    
      CkArrayOptions opts(numBlks, numBlks);
      switch (mapping) {
      case 1:
        opts.setMap(CProxy_LUBalancedSnakeMap2::ckNew(numBlks, BLKSIZE));
	break;
      case 0:
	opts.setMap(CProxy_BlockCyclicMap::ckNew());
	break;
      }
      luArrProxy = CProxy_LUBlk::ckNew(opts);

      workStarted = true;
    
      luArrProxy.init(0, BLKSIZE, numBlks, memThreshold);
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

  void terminateProg(CkReductionMsg *msg) {
    CkCallback cb(CkIndex_Main::done(NULL),thisProxy); 
    traceCriticalPathBack(cb);
  }

  void done(pathInformationMsg *m){
    // CkPrintf("Main::done() After critical path has been determined\n");
    //	  m->printme();
    gotoNextPhase(); // <<< Make sure we get timings for the phase that just finished.
    CkExit();
  }

};

class LUBlk: public CBase_LUBlk {
  /// Variables used during factorization
  double *LU;
  //VALIDATION: variable to hold copy of untouched matrix
  double *A;
  int BLKSIZE, numBlks;
  blkMsg *L, *U;
  int internalStep, activeCol;

  /// Variables used only during solution
  double *bvec;
  //VALIDATION: variable to hold copy of untouched b vector
  double *b;
  //VALIDATION: variable to hold copy of final x
  double *x;

  double* storedVec;
  int diagRec;

  //VALIDATION: count of the number of point-to-point messages rcvd in a row
  int msgsRecvd;

  //Variables for pivoting SDAG code

  int row1Index, row2Index, localRow1, localRow2,
    otherRowIndex, thisLocalRow, globalThisRow, globalOtherRow;
  bool remoteSwap;
  locval l;
  int pivotBlk;

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
	  // variable A has the original matrix sub-block
	  // solution sub-vector x is in variable x on the diagonals
	  // variable b has the original b vector on the diagonals

	  //Diagonals distribute x across entire column
	  if(thisIndex.x == thisIndex.y)
	  {
		  CProxySection_LUBlk col =
				  CProxySection_LUBlk::ckNew(thisArrayID, 0, numBlks-1,
						  1, thisIndex.y, thisIndex.y, 1);

		  col.recvXvec(BLKSIZE, x);
	  }
  }

  //VALIDATION
  void recvXvec(int size, double* xvec) {

	  double *partial_b = new double[BLKSIZE];

	  //Perform local dgemv
#if USE_ESSL
    dgemv("T", BLKSIZE, BLKSIZE, 1.0, A, BLKSIZE, x, 1, 0.0, partial_b, 1);
#else
    cblas_dgemv( CblasColMajor, CblasTrans,
    		  BLKSIZE, BLKSIZE, 1.0, A,
    		  BLKSIZE, xvec, 1, 0.0, partial_b, 1);
#endif

    //sum-reduction of result across row with diagonal element as target
    thisProxy(thisIndex.x,thisIndex.x).sumBvec(BLKSIZE,partial_b);

    //if you are not the diagonal, signal finish
    if(thisIndex.x != thisIndex.y)
    	contribute(CkCallback(CkIndex_Main::terminateProg(NULL),mainProxy));
  }

  //VALIDATION
  void sumBvec(int size, double* partial_b) {

	  //Clear bvec before first message processed for sum-reduction
	  if(msgsRecvd == 0) {
		  for(int i = 0; i < BLKSIZE; i++)
			  bvec[i] = 0.0;
	  }

	  //Sum up messages
	  if(++msgsRecvd <= numBlks) {
		  for (int i = 0; i < size; i++) {
			  bvec[i] += partial_b[i];
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
		  residuals[i] = b[i] - bvec[i];
		  if(fabs(residuals[i]) > 1e-10)
			  CkPrintf("WARNING: RESIDUAL VALUE = %e\n",residuals[i]);
		  //CkPrintf("res[%d] = %f - %f = %e\n",i,b[i],bvec[i],residuals[i]);
	  }

	  //Diagonal elements signal finish
	  contribute(CkCallback(CkIndex_Main::terminateProg(NULL),mainProxy));
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
      for (int i=0; i<blocksize*blocksize; i++) { 
	m1[i] = rnd.toRndDouble(rnd.nextRndInt()); 
	m2[i] = rnd.toRndDouble(rnd.nextRndInt()); 
	m3[i] = rnd.toRndDouble(rnd.nextRndInt()); 
      } 
      
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

  void initVec(int size, double* vec) {
    bvec = new double[BLKSIZE];
    memcpy(bvec, vec, sizeof(double) * BLKSIZE);

    //VALIDATION: keep a copy of the original b vector on the diagonals
    if(thisIndex.x == thisIndex.y)
    {
    	b = new double[BLKSIZE];
    	memcpy(b, bvec, sizeof(double) * BLKSIZE);
    }

#if defined(PRINT_VECTORS)
    for (int i = 0; i < BLKSIZE; i++) {
      CkPrintf("memcpy bvec[%d] = %f\n", i, bvec[i]);
    }
#endif

    contribute(CkCallback(CkIndex_Main::finishInit(), mainProxy));
  }

  void init(int _whichMulticastStrategy, int _BLKSIZE, int _numBlks, int memThreshold) {
    whichMulticastStrategy = _whichMulticastStrategy;
    BLKSIZE = _BLKSIZE;
    numBlks = _numBlks;
    
    // Set the schedulers memory usage threshold to the one based upon a control point
    schedAdaptMemThresholdMB = memThreshold;

    CkAssert(BLKSIZE>0); // If this fails, readonly variables aren't
			 // propagated soon enough. I'm assuming they
			 // are safe to use here.

#if USE_MEMALIGN
    LU = (double*)memalign(128, BLKSIZE*BLKSIZE*sizeof(double) );
    //VALIDATION: A is an untouched copy for validation
    A  = (double*)memalign(128, BLKSIZE*BLKSIZE*sizeof(double) );
    //	 CkPrintf("LU mod 128 = %lu\n", ((unsigned long)LU) % 128);
    CkAssert(LU != NULL);
#else
    LU = new double[BLKSIZE*BLKSIZE];
    //VALIDATION: A is an untouched copy for validation
    A  = new double[BLKSIZE*BLKSIZE];
#endif

    internalStep = 0;  
     
    traceUserSuppliedData(-1);	
    traceMemoryUsage();	 
     
    //MatGen rnd(thisIndex.x * numBlks + thisIndex.y);

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

    //VALIDATION: Make a copy for validation
    memcpy(A, LU, BLKSIZE*BLKSIZE*sizeof(double));

    this->print("input-generated-LU");

    testdgemm();

    contribute(CkCallback(CkIndex_Main::finishInit(), mainProxy));
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


  //thisIndex.x indicates the internal step it is going to work on
  void solveLocalLU() {
    traceLU t(internalStep, traceSolveLocalLU);

    CkAssert(internalStep == thisIndex.x);

    DEBUG_PRINT("elem[%d,%d]::solveLocalLU called at step %d\n", thisIndex.x, thisIndex.y, internalStep);

    CkAssert(thisIndex.x == thisIndex.y);

    LUdecompose(LU);
  }


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
    dtrsm("L", "L", "N", "U", BLKSIZE, BLKSIZE, 1.0, givenL, BLKSIZE, LU, BLKSIZE);
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
	   -1.0, incomingL,
	   BLKSIZE, incomingU, BLKSIZE,
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

  void processComputeL(int ignoredParam) {
    DEBUG_PRINT("processComputeL() called on block %d,%d\n", thisIndex.x, thisIndex.y);
    CkAssert(internalStep==thisIndex.y && U);
    
    // We are in the left row of active blocks
      
    //  CkPrintf("[%d] chare %d,%d internalStep=%d computeL\n", CkMyPe(), thisIndex.x, thisIndex.y, internalStep);
    computeL(U);
      
    //broadcast the newly computed L rightwards to the blocks in the same row
    DEBUG_PRINT("[%d] chare %d,%d is left block this step, multicast L\n", CkMyPe(), thisIndex.x, thisIndex.y);
    multicastRecvL();
      
    CmiFree(UsrToEnv(U));

    DEBUG_PRINT("chare %d,%d is now done\n",  thisIndex.x, thisIndex.y);
  }
  
  void processLocalLU(int ignoredParam) {
    DEBUG_PRINT("processLocalLU() called on block %d,%d\n", thisIndex.x, thisIndex.y);
    // We are the top-left-most active block
    CkAssert(internalStep==thisIndex.x && internalStep==thisIndex.y);

//    double mem = CmiMemoryUsage();
//    CkPrintf("CmiMemoryUsage() = %lf\n", mem);

    // CkPrintf("[%d] chare %d,%d internalStep=%d solveLocalLU\n", CkMyPe(), thisIndex.x, thisIndex.y, internalStep);
    solveLocalLU(); // compute a local LU

    // If this is the very last bottom rightmost block
    if (thisIndex.x == numBlks-1 && thisIndex.y == numBlks-1) {
      DEBUG_PRINT("[%d] chare %d,%d calling mainProxy.iterationCompleted\n", CkMyPe(), thisIndex.x, thisIndex.y);
      mainProxy.iterationCompleted();
    } else {
      DEBUG_PRINT("[%d] chare %d,%d is top,left block this step, multicast U & L\n", CkMyPe(), thisIndex.x, thisIndex.y);
      multicastRecvU();	//broadcast the U downwards to the blocks in the same column
      multicastRecvL();		//broadcast the L rightwards to the blocks in the same row
    }
      
    DEBUG_PRINT("chare %d,%d is now done\n",  thisIndex.x, thisIndex.y);
  }

#if 0
  /// Call progress on myself, possibly using priorities 
  inline void selfContinue(){
    int integerPrio;
    
    //	  CkPrintf("continuing %d,%d  internalStep=%d \n", thisIndex.x,thisIndex.y, internalStep);


    double c1 = 1.0;
    double c2 = 1.0;
    double c3 = 1.0;
    double c4 = 1.0;

#if 1
    // Low priority trailing updates
    // High priorities for critical path (solve local LUs)
    if(thisIndex.x == thisIndex.y && thisIndex.x == internalStep){
      integerPrio = -1000; // highest priority
    } else if(thisIndex.x == internalStep || thisIndex.y == internalStep){
      integerPrio = c1*(-1*(BLKSIZE-internalStep) -5) + c2;
    } else {
      // Trailing updates have lower priorities that increase from top left to bottom right
      integerPrio = (internalStep+1)*c3 + (thisIndex.x+thisIndex.y)*c4;
    }
#elif 0
    // Low priority trailing updates
    // High priorities for critical path (solve local LUs)
    if(thisIndex.x == thisIndex.y && thisIndex.x == internalStep){
      integerPrio = -10*c1; // highest priority
    } else if(thisIndex.x == internalStep || thisIndex.y == internalStep){
      integerPrio = -5; 
    } else if(thisIndex.x == thisIndex.y){
      integerPrio = -1; // high priority
    } else {
      // Trailing updates have lower priorities that increase from top left to bottom right
      integerPrio = (internalStep+1)*c3 + (thisIndex.x+thisIndex.y)*c4;
    }
#elif 0
    // High priorities for trailing updates
    if(thisIndex.x == thisIndex.y && thisIndex.x == internalStep){
      integerPrio = 10; // corners
    } else if(thisIndex.x == internalStep || thisIndex.y == internalStep){
      integerPrio = 9; // edges
    } else {
      // Trailing updates
      integerPrio = -100;
    }
#else
    // High priorities for early trailing updates
    if(thisIndex.x == thisIndex.y && thisIndex.x == internalStep){
      integerPrio = -10000; // highest priority
    } else if(thisIndex.x == internalStep || thisIndex.y == internalStep){
      integerPrio = c1*(internalStep-BLKSIZE -5) + c2;
    } else {
      // Trailing updates have lower priorities that increase from top left to bottom right
      integerPrio = (internalStep+1)*c3 + (thisIndex.x+thisIndex.y)*c4;
      if(internalStep < 5){
	  integerPrio -= (5-internalStep)*(BLKSIZE*3);
      }

    }

#endif	  

    CkEntryOptions eOpts; 
    eOpts.setPriority (integerPrio); // setPriority sets the queuing type internally
    
    
    switch(canContinue()) {
    case NO_CONTINUE:
      break;
    case CONTINUE_LU:
      thisProxy(thisIndex.x,thisIndex.y).processLocalLU(0, &eOpts);
      break;
    case CONTINUE_U:
      thisProxy(thisIndex.x,thisIndex.y).processComputeU(0, &eOpts);
      break;
    case CONTINUE_L:
      thisProxy(thisIndex.x,thisIndex.y).processComputeL(0, &eOpts);
      break;
    case CONTINUE_TRAIL:
      thisProxy(thisIndex.x,thisIndex.y).processTrailingUpdate(0, &eOpts);
      break;  
    }

  }
#endif

  int getIndex(int i, int j) {
    return i * BLKSIZE + j;
  }

  void localSolve(double *xvec, double *preVec, bool diag, bool forward) {
    if (!diag) {
      for (int i = 0; i < BLKSIZE; i++) {
        xvec[i] = 0.0;
      }
      
      for (int i = 0; i < BLKSIZE; i++) {
        for (int j = 0; j < BLKSIZE; j++) {
          xvec[i] += LU[getIndex(i,j)] * preVec[j];
        }
      }
    } else {
      if (preVec != NULL)
        for (int i = 0; i < BLKSIZE; i++)
          xvec[i] = bvec[i] - preVec[i];
      else
        memcpy(xvec, bvec, sizeof(double) * BLKSIZE);
      
      if (forward)
	localForward(xvec);
      else
	localBackward(xvec);
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

  void solve(bool backward, int size, double* preVec) {
    CkPrintf("solved called on: (%d, %d)\n", thisIndex.x, thisIndex.y);

    bool forward = !backward;

#if defined(PRINT_VECTORS)
    for (int i = 0; i < BLKSIZE; i++) {
      CkPrintf("bvec[%d] = %f\n", i, bvec[i]);
    }
#endif

    CkPrintf("allocate x\n");
    x = new double[BLKSIZE];

    CkPrintf("Calling local solver - diag = true\n");
    localSolve(x, (size == BLKSIZE) ? preVec : NULL, true, forward);

    if (forward) {
      
      /*for (int i = 0; i < BLKSIZE; i++) {
        CkPrintf("%d, %d: xvec[%d] = %f\n", thisIndex.x, thisIndex.y, i, xvec[i]);
        }*/

      for (int i = 0; i < BLKSIZE; i++) {
        bvec[i] = x[i];
#if defined(PRINT_VECTORS)
        CkPrintf("After forward solve - x[%d] = %f\n", i, x[i]);
#endif
      }

      delete[] storedVec;
      storedVec = NULL;
      
      diagRec = 0;
      
      if (thisIndex.x == numBlks-1) {
	CkPrintf("forward-solve complete, (%d, %d) numBlks = %d\n", thisIndex.x, thisIndex.y, numBlks);

        // Perform backward solving
        thisProxy(numBlks-1, numBlks-1).solve(true, 0, NULL);
      }

      CProxySection_LUBlk col = 
	CProxySection_LUBlk::ckNew(thisArrayID, thisIndex.x+1, numBlks-1, 
				   1, thisIndex.y, thisIndex.y, 1);

      col.forwardSolve(BLKSIZE, x);
    } else {
#if defined(PRINT_VECTORS)
      for (int i = 0; i < BLKSIZE; i++) {
        CkPrintf("After backward solve - x[%d] = %f\n", i, x[i]);
      }
#endif

      //VALIDATION: put final solution into bvec
//      for (int i = 0; i < BLKSIZE; i++)
//      	bvec[i] = x[i];

      if (thisIndex.x == 0) {
    	  CkPrintf("backward-solve complete, (%d, %d) numBlks = %d\n", thisIndex.x, thisIndex.y, numBlks);

    	  mainProxy.iterationCompleted();
      } else {
    	  CProxySection_LUBlk col =
    			  CProxySection_LUBlk::ckNew(thisArrayID, 0, thisIndex.x-1,
    					  1, thisIndex.y, thisIndex.y, 1);

    	  col.backwardSolve(BLKSIZE, x);
      }
    }
  }

  void diagForwardSolve(int size, double* vec) {
    if (!storedVec) {
      storedVec = new double[BLKSIZE];
      for (int i = 0; i < BLKSIZE; i++) {
	storedVec[i] = 0.0;
      }
    }

    for (int i = 0; i < BLKSIZE; i++) {
      storedVec[i] += vec[i];
    }
    
    CkPrintf("FORWARD: expected %d, received %d\n", thisIndex.y, diagRec+1);

    if (++diagRec == thisIndex.y)
      thisProxy(thisIndex.x, thisIndex.y).solve(false, BLKSIZE, storedVec);
  }

  void diagBackwardSolve(int size, double* vec) {
    if (!storedVec) {
      storedVec = new double[BLKSIZE];
      for (int i = 0; i < BLKSIZE; i++) {
	storedVec[i] = 0.0;
      }
    }

    for (int i = 0; i < BLKSIZE; i++) {
      storedVec[i] += vec[i];
    }
    
    int expected = numBlks-1-thisIndex.y;
    
    CkPrintf("BACKWARD: expected %d, received %d\n", expected, diagRec+1);

    if (++diagRec == expected)
      thisProxy(thisIndex.x, thisIndex.y).solve(true, BLKSIZE, storedVec);
  }

  void forwardSolve(int size, double* preVec) {
    double *xvec = new double[BLKSIZE];
    
    localSolve(xvec, preVec, false, true);
    
    CkPrintf("diagForwardSolve called from: (%d, %d), on: (%d, %d)\n", 
             thisIndex.x, thisIndex.y, thisIndex.x, thisIndex.x);
    
    thisProxy(thisIndex.x, thisIndex.x).diagForwardSolve(size, xvec);
  }

  void backwardSolve(int size, double* preVec) {
    double *xvec = new double[BLKSIZE];
    
    localSolve(xvec, preVec, false, false);
    
    CkPrintf("diagBackwardSolve called from: (%d, %d), on: (%d, %d)\n", 
             thisIndex.x, thisIndex.y, thisIndex.x, thisIndex.x);
    
    thisProxy(thisIndex.x, thisIndex.x).diagBackwardSolve(size, xvec);
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

  void LUdecompose(double* A) { 
    for (int j = 0; j < BLKSIZE; j++) {
      for (int i = 0; i <= j; i++) {
	double sum = 0.0;
	for (int k = 0; k < i; k++)
	  sum += A[getIndex(i, k)] * A[getIndex(k,j)];
	A[getIndex(i,j)] -= sum;
      }
      for (int i = j + 1; i < BLKSIZE; i++) {
	double sum = 0.0;
	for (int k = 0; k < j; k++)
	  sum += A[getIndex(i,k)] * A[getIndex(k,j)];
	A[getIndex(i,j)] -= sum;
	A[getIndex(i,j)] /= A[getIndex(j,j)];
      }
    }
  }

  //Performs data exchange
  void applySwap(int loc, double* data) {

  }

  //Should exchange local data (called by applySwap?)
  void swapLocal(int loc, int col) {

  }

  //Does the local multiplier computation after U is sent to the blocks below??
  void diagonalUpdate(int col) {
	  computeMultipliers(LU[getIndex(col,col)],col,col);

	  for(int k=col+1;k<BLKSIZE;k++) {
		  for(int j=k; j<BLKSIZE; j++) {
			  //U[k] might need to be U[j]?
			  LU[getIndex(j,k)] =  LU[getIndex(j,k)] - LU[getIndex(j,col)]*U[k];
		  }
	  }
  }

private:
  //internal functions for creating messages to encapsulate the priority
  inline blkMsg* createABlkMsg() {
    blkMsg *msg;
    
    if(doPrioritize) {
      msg = new(BLKSIZE*BLKSIZE, 8*sizeof(int)) blkMsg;
      DEBUG_PRINT("setting priority to internalStep=%d\n", internalStep);
      *((int*)CkPriorityPtr(msg)) = (int)internalStep;
      CkSetQueueing(msg, CK_QUEUEING_IFIFO);
    } else {
      msg = new(BLKSIZE*BLKSIZE) blkMsg;
    }
    msg->setMsgData(LU, internalStep, BLKSIZE);

    return msg;
  }

  locval findLocVal(int startRow, int col, double curMax, int curRowMax) {
    locval l;
    l.val = curMax;
    l.loc = curRowMax;
    for (int row = startRow; row < BLKSIZE; row++) {
      if (LU[getIndex(row, col)] > l.val) {
	l.val = LU[getIndex(row, col)];
	l.loc = row + BLKSIZE * thisIndex.y;
      }
    }
    return l;
  }


  /// Compute the multipliers based on the pivot value from the diagonal chare
  //  starting at [row, col]
  void computeMultipliers(double a_kk, int row, int col) {
	  for(int i = row; i<BLKSIZE;i++)
		  LU[getIndex(i,col)] = LU[getIndex(i,col)]/a_kk;
  }

  /// Update the values in the columns ahead of the active column within the
  // pivot section based on the Usegment and the multipliers
  void updateAllCols(int col, double* U) {

	  //TODO: Replace with DGEMM
	  for(int k=col+1;k<BLKSIZE;k++) {
		  for(int j=0; j<BLKSIZE; j++) {
			  //U[k] might need to be U[j]?
			  LU[getIndex(j,k)] =  LU[getIndex(j,k)] - LU[getIndex(j,col)]*U[k];
		  }
	  }
  }

};

#include "lu.def.h"
