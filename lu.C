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
//#include <string.h>
#include <iostream>
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


/// Called periodically for debugging
static void periodicDebug(void* ptr, double currWallTime){
  //  CcdCallFnAfterOnPE((CcdVoidFn)periodicDebug, (void*)NULL, controlPointSamplePeriod, CkMyPe());
  luArrProxy.printInfo();
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

class blkMsg: public CMessage_blkMsg {
public:
  int step;
  char pad[16-((sizeof(envelope)+sizeof(int))%16)];
  double *data;
  int BLKSIZE;

  blkMsg(int _BLKSIZE) : BLKSIZE(_BLKSIZE) { }

  void setMsgData(double *d, int s) {
    memcpy(data, d, sizeof(double)*BLKSIZE*BLKSIZE);
    step = s;
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

  Main(CkArgMsg* m) : numIterations(1) {
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
 
    if (gMatSize%1024!=0) 
      CkAbort("The matrix size should be a multiple of 1024!\n");
  
    CkPrintf("Running LU on %d processors (%d nodes) on matrix %dX%d with control points\n",
	     CkNumPes(), CmiNumNodes(), gMatSize, gMatSize);

    multicastStats[0] = ComlibRegister(new OneTimeRingMulticastStrategy() ); 
    multicastStats[1] = ComlibRegister(new OneTimeNodeTreeMulticastStrategy(2) ); 
    multicastStats[2] = ComlibRegister(new OneTimeNodeTreeMulticastStrategy(3) ); 
    multicastStats[3] = ComlibRegister(new OneTimeNodeTreeMulticastStrategy(4) ); 

    ControlPoint::EffectIncrease::GrainSize("block_size");
    ControlPoint::EffectDecrease::Concurrency("block_size");
    
    //    ControlPoint::EffectIncrease::Concurrency("mapping");
    ControlPoint::EffectIncrease::NumMessages("mapping");
    ControlPoint::EffectIncrease::MessageOverhead("mapping");

    ControlPoint::EffectIncrease::MemoryConsumption("memory_threshold");


    thisProxy.iterationCompleted();

  }

  void finishInit() {
    //    CkPrintf("finishInit called\n");
    luArrProxy.flushLogs();
  }

  void continueIter() {
    CkCallback *cb = new CkCallback(CkIndex_Main::arrayIsCreated(NULL), thisProxy);
    CkStartQD(*cb); // required currently for use with new Comlib
  }

  void arrayIsCreated(CkReductionMsg * m) {
    startTime = CmiWallTimer();
    luArrProxy(0,0).processLocalLU(0);
  }
  


  void iterationCompleted() {
    if(iteration == 0){
      CkPrintf("Initialization Complete\n");
    } else {
      double endTime = CmiWallTimer();
      double duration = endTime - startTime;
      registerControlPointTiming(duration);
      
      CkPrintf("Iteration %d time: %fs\n", iteration, duration);
      
      if (iteration == numIterations){ 

	outputStats();

	terminateProg();
	return;
      }
      
      luArrProxy.ckDestroy();
    }


    iteration++;
    
    // Only advance phases after a few factorizations have been performed
    // Prior to the first phase of actual work, iteration=1
    if( 1 || iteration % 2 == 1 || iteration==1){
      gotoNextPhase();
      
      whichMulticastStrategy = controlPoint("multicast_strategy", 2, 2);
      BLKSIZE = 1 << controlPoint("block_size", 9, 9);
      mapping = controlPoint("mapping", 1, 1);
      memThreshold = 200 + controlPoint("memory_threshold", 0, 12) * 100;
      
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
    
    switch (mapping) {
    case 1: {
      CProxy_LUBalancedSnakeMap2 map = CProxy_LUBalancedSnakeMap2::ckNew(numBlks, BLKSIZE);
      CkArrayOptions opts(numBlks, numBlks);
      opts.setMap(map);
      luArrProxy = CProxy_LUBlk::ckNew(opts);
    } break;
    case 0: {
      CProxy_BlockCyclicMap map = CProxy_BlockCyclicMap::ckNew();
      CkArrayOptions opts(numBlks, numBlks);
      opts.setMap(map);
      luArrProxy = CProxy_LUBlk::ckNew(opts);
    } break;
    }
    
    luArrProxy.init(0, BLKSIZE, numBlks, memThreshold);
  }
  
  void outputStats() {
    double endTime = CmiWallTimer();
    double duration = endTime-startTime;

    double n = gMatSize;

    long long flopCount = 0;     // floating point ops
    for (int i=1;i<=gMatSize;i++) {
      for (int j=1+i; j<=gMatSize; j++) {
	flopCount += (1+2*gMatSize-2*i);
      }
    }

    double flops = ((double)flopCount)  / duration; // floating point ops per second
    double gflops = flops / 1000000000.0; // Giga fp ops per second
    std::cout << "RESULT procs: \t" << CkNumPes() << "\tblock size:\t" << BLKSIZE << "\tGFlops:\t" << gflops << "\tTime(s):\t" << duration << std::endl;

    double HPL_flop_count =  (2.0/3.0*n*n*n-2*n*n)/duration ;
    double HPL_gflops =  HPL_flop_count / 1000000000.0; // Giga fp ops per second
    std::cout << "HPL flop count gives \t" << HPL_gflops << "\tGFlops" << std::endl;

    
    double gflops_per_core = HPL_gflops / (double)CkNumPes();
    double fractionOfPeakOnOrder = gflops_per_core / 7.4585;
    double fractionOfPeakOnOrderPercent = fractionOfPeakOnOrder * 100.0;
    
    double fractionOfPeakOnAbe = gflops_per_core / 9.332;
    double fractionOfPeakOnAbePercent = fractionOfPeakOnAbe * 100.0;

    double fractionOfPeakOnKraken =  gflops_per_core / 10.4;
    double fractionOfPeakOnKrakenPercent = fractionOfPeakOnKraken * 100.0;
    
    double fractionOfPeakOnBGP =  gflops_per_core / 3.4;
    double fractionOfPeakOnBGPPercent = fractionOfPeakOnBGP * 100.0;

    std::cout << "If ran on order.cs.uiuc.edu, I think you got  \t" << fractionOfPeakOnOrderPercent << "% of peak" << std::endl;
    std::cout << "If ran on abe.ncsa.uiuc.edu, I think you got  \t" << fractionOfPeakOnAbePercent << "% of peak" << std::endl;
    std::cout << "If ran on kraken, I think you got  \t" << fractionOfPeakOnKrakenPercent << "% of peak" << std::endl;
    std::cout << "If ran on BG/P, I think you got  \t" << fractionOfPeakOnBGPPercent << "% of peak" << std::endl;

  }

  void terminateProg() {
    CkCallback cb(CkIndex_Main::done(NULL),thisProxy); 
    traceCriticalPathBack(cb);
  }

  void done(pathInformationMsg *m){
    // CkPrintf("Main::done() After critical path has been determined\n");
    //    m->printme();
    gotoNextPhase(); // <<< Make sure we get timings for the phase that just finished.
    CkExit();
  }

};

class LUBlk: public CBase_LUBlk {
public:
  double *LU;
  int whichMulticastStrategy;
  bool done;
  int alreadyReEnqueuedDuringPhase;
  int BLKSIZE, numBlks;
  
private:
  CkVec<blkMsg *> UBuffers;
  CkVec<blkMsg *> LBuffers;

  int internalStep;


  MERGE_PATH_DECLARE_D(A);


public:
  LUBlk() {

    whichMulticastStrategy = 0;
    done = false;
    alreadyReEnqueuedDuringPhase = -1;

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

  void init(int _whichMulticastStrategy, int _BLKSIZE, int _numBlks, int memThreshold){
    whichMulticastStrategy = _whichMulticastStrategy;
    BLKSIZE = _BLKSIZE;
    numBlks = _numBlks;
    
    // Set the schedulers memory usage threshold to the one based upon a control point
    schedAdaptMemThresholdMB = memThreshold;

    done = false;
    alreadyReEnqueuedDuringPhase = -1;

    CkAssert(BLKSIZE>0); // If this fails, readonly variables aren't
			 // propagated soon enough. I'm assuming they
			 // are safe to use here.

#if USE_MEMALIGN
    LU = (double*)memalign(128, BLKSIZE*BLKSIZE*sizeof(double) );
    //   CkPrintf("LU mod 128 = %lu\n", ((unsigned long)LU) % 128);
    CkAssert(LU != NULL);
#else
    LU = new double[BLKSIZE*BLKSIZE];
#endif

    internalStep = 0;  
     
    traceUserSuppliedData(-1);  
    traceMemoryUsage();  
     
    MatGen rnd(thisIndex.x * numBlks + thisIndex.y);      
    for (int i=0; i<BLKSIZE*BLKSIZE; i++) {  
      LU[i] = rnd.toRndDouble(rnd.nextRndInt());  
    } 


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
//     CkAssert(!done);
//     thisProxy(0,0).processLocalLU(0);
//   }


  /* Computation functions that should be called localy related with each
   * state of the block
   */


  //thisIndex.x indicates the internal step it is going to work on
  void solveLocalLU() {
    traceUserSuppliedData(internalStep);
    traceMemoryUsage();
    internalStep = thisIndex.x;

    DEBUG_PRINT("elem[%d,%d]::solveLocalLU called at step %d\n", thisIndex.x, thisIndex.y, internalStep);

    CkAssert(thisIndex.x == thisIndex.y);

    double luStart = CmiWallTimer();

    //#ifdef USE_LAPACK

    //There's an output for permuation array which is not
    //used in the current non-lapack version. This permuation
    //array should also be broadcasted to those elements that
    //only needs ComputeL and ComputeU (???)

#if USE_MKL_CBLAS_H 
    int size = BLKSIZE;
    int *ipiv = new int[BLKSIZE];
    int info;
    // This one doesn't quite do what we want... it does pivoting
    dgetrf(&size, &size, LU, &size, ipiv, &info);
#elif USE_ACCELERATE_BLAS
    __CLPK_integer size = static_cast<__CLPK_integer>(BLKSIZE);
    __CLPK_integer *ipiv = new __CLPK_integer[BLKSIZE];
    __CLPK_integer info;
    dgetrf_(&size, &size, LU, &size, ipiv, &info);
    delete[] ipiv;
#elif USE_ESSL
    int size = BLKSIZE;
    int *ipiv = new int[BLKSIZE];
    int info;
    // This one doesn't quite do what we want... it does pivoting
    dgetrf(size, size, LU, size, ipiv, &info);
    delete [] ipiv;
#else
    int *ipiv = new int[BLKSIZE];
    clapack_dgetrf(CblasRowMajor, BLKSIZE, BLKSIZE, LU, BLKSIZE, ipiv);


    //    int clapack_dgetrf(const enum CBLAS_ORDER Order, const int
    //M, const int N,
    //	       double *A, const int lda, int *ipiv);    

    delete [] ipiv;
#endif
    
    /** @FIXME: do the permutation of the rows specified by ipiv */

    traceUserBracketEvent(traceSolveLocalLU, luStart, CmiWallTimer());
  }


  void computeU(blkMsg *givenLMsg) {
    traceUserSuppliedData(internalStep);
    traceMemoryUsage();
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

    double uStart = CmiWallTimer();

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

    traceUserBracketEvent(traceComputeU, uStart, CmiWallTimer());
  }

  void computeL(blkMsg *givenUMsg) {
    traceUserSuppliedData(internalStep);
    traceMemoryUsage();
    double *givenU = givenUMsg->data;
    //    CkAssert( ((unsigned long)givenU) % 16 == 0);


    //traceUserEvent(traceComputeL);
    double lStart = CmiWallTimer();

    DEBUG_PRINT("elem[%d,%d]::computeL called at step %d\n", thisIndex.x, thisIndex.y, internalStep);

#if USE_ESSL
    dtrsm("R", "U", "N", "N", BLKSIZE, BLKSIZE, 1.0, givenU, BLKSIZE, LU, BLKSIZE);
#else
    cblas_dtrsm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, BLKSIZE, BLKSIZE, 1.0, givenU, BLKSIZE, LU, BLKSIZE);
#endif

    traceUserBracketEvent(traceComputeL, lStart, CmiWallTimer());
  }

  void updateMatrix(blkMsg *givenLMsg, blkMsg *givenUMsg) {
    traceUserSuppliedData(internalStep);
    traceMemoryUsage();
    DEBUG_PRINT("elem[%d,%d] is updating its value at step %d\n", thisIndex.x, thisIndex.y, internalStep);

    //traceUserEvent(traceTrailingUpdate);
    double updateStart = CmiWallTimer();

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

    traceUserBracketEvent(traceTrailingUpdate, updateStart, CmiWallTimer());
  }


  void updateRecvU(blkMsg *UMsg) {
    CkAssert(!done);
    //   CkAssert(UMsg->step == 0);
    traceUserSuppliedData(UMsg->step);
    traceMemoryUsage();
    DEBUG_PRINT("elem[%d,%d]::updateRecvU entry method containing message from step %d\n", thisIndex.x, thisIndex.y, UMsg->step);
    bufferU(UMsg);
    MERGE_PATH_MAX_D(A,UMsg->step);
    if(canContinue() && alreadyReEnqueuedDuringPhase < internalStep){
      alreadyReEnqueuedDuringPhase = internalStep;
      DEBUG_PRINT("NOTE   :                                  calling progress() from updateRecvU\n");
      selfContinue();
    }
  }

  void updateRecvL(blkMsg *LMsg) {
    CkAssert(!done);
    traceUserSuppliedData(LMsg->step);
    traceMemoryUsage();
    DEBUG_PRINT("elem[%d,%d]::updateRecvL entry method containing message from step %d\n", thisIndex.x, thisIndex.y, LMsg->step);
    bufferL(LMsg);
    MERGE_PATH_MAX_D(A,LMsg->step);
    if(canContinue() && alreadyReEnqueuedDuringPhase < internalStep){
      alreadyReEnqueuedDuringPhase = internalStep;
      DEBUG_PRINT("NOTE   :                                  calling progress() from updateRecvL\n");
      selfContinue();
    }
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
    givenU->setMsgData(LU, internalStep);
    oneCol.updateRecvU(givenU);

//     for(int i=thisIndex.x+1; i<numBlks; i++){
//       blkMsg *givenU = createABlkMsg();
//       givenU->setMsgData(LU, internalStep);
//       DEBUG_PRINT("P2P sending U from %d,%d down to %d,%d\n", thisIndex.x, thisIndex.y, i,thisIndex.y);
//       thisProxy(i,thisIndex.y).updateRecvU(givenU);
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
    givenL->setMsgData(LU, internalStep);
    //CkAssert(givenL->step == 0);
    oneRow.updateRecvL(givenL);
    
//     for(int i=thisIndex.y+1; i<numBlks; i++){
//       blkMsg *givenL = createABlkMsg();
//       givenL->setMsgData(LU, internalStep);
//       DEBUG_PRINT("P2P sending L from %d,%d right to %d,%d\n", thisIndex.x, thisIndex.y, thisIndex.x, i);
//       thisProxy(thisIndex.x, i).updateRecvL(givenL);
//     }
    
  }

  void processComputeU(int ignoredParam) {
    DEBUG_PRINT("processComputeU() called on block %d,%d\n", thisIndex.x, thisIndex.y);
    CkAssert(!done);
    CkAssert(internalStep==thisIndex.x && getBufferedL(internalStep)!=NULL);
    // We are in the top row of active blocks, and we
    // have received the incoming L
	
    // CkPrintf("[%d] chare %d,%d internalStep=%d computeU\n", CkMyPe(), thisIndex.x, thisIndex.y, internalStep);
    computeU(LBuffers[internalStep]);
    
    DEBUG_PRINT("[%d] chare %d,%d is top block this step, multicast U\n", CkMyPe(), thisIndex.x, thisIndex.y);
    multicastRecvU(); //broadcast the newly computed U downwards to the blocks in the same column
    
    removeBufferedL(internalStep);	// Cleanup
    
    // Verify that there are no outstanding buffered messages.
    CkAssert(buffersEmpty());
    
    deallocateBuffers();
    MERGE_PATH_DELETE_D(A,internalStep);
    // This block is now done
    DEBUG_PRINT("chare %d,%d is now done\n",  thisIndex.x, thisIndex.y);
    done = true;
    
  }
  
  void processComputeL(int ignoredParam) {
    DEBUG_PRINT("processComputeL() called on block %d,%d\n", thisIndex.x, thisIndex.y);
    CkAssert(!done);
    CkAssert(internalStep==thisIndex.y && getBufferedU(internalStep)!=NULL);
    
    // We are in the left row of active blocks
      
    //  CkPrintf("[%d] chare %d,%d internalStep=%d computeL\n", CkMyPe(), thisIndex.x, thisIndex.y, internalStep);
    computeL(UBuffers[internalStep]);
      
    //broadcast the newly computed L rightwards to the blocks in the same row
    DEBUG_PRINT("[%d] chare %d,%d is left block this step, multicast L\n", CkMyPe(), thisIndex.x, thisIndex.y);
    multicastRecvL();
      
    removeBufferedU(internalStep);	// Cleanup

    // Verify that there are no outstanding buffered messages.
    CkAssert(buffersEmpty());

    deallocateBuffers();
    MERGE_PATH_DELETE_D(A,internalStep);
    // This block is now done
    DEBUG_PRINT("chare %d,%d is now done\n",  thisIndex.x, thisIndex.y);
    done = true;
    
  }
  
  void processTrailingUpdate(int ignoredParam) {
    DEBUG_PRINT("processTrailingUpdate() called on block %d,%d\n", thisIndex.x, thisIndex.y);
    CkAssert(!done);
    
    double *incomingL = getBufferedL(internalStep);
    double *incomingU = getBufferedU(internalStep);

    CkAssert(incomingL!=NULL && incomingU!=NULL);
    
    // Do trailing update
    
    //CkPrintf("[%d] chare %d,%d internalStep=%d updateMatrix\n", CkMyPe(), thisIndex.x, thisIndex.y, internalStep);
    updateMatrix(LBuffers[internalStep], UBuffers[internalStep]);
    
    removeBufferedL(internalStep);
    removeBufferedU(internalStep);
    
    DEBUG_PRINT("elem[%d,%d] advancing to step %d\n", thisIndex.x, thisIndex.y, internalStep+1);
    internalStep++;
    

    // If we have received the L and U messages out of order, then some of the 
    // calls into this progress method will not do anything useful.
    // Thus we need to tell ourself to continue if possible, because
    // we cannot rely upon a future progress call being made.
    if(canContinue() && alreadyReEnqueuedDuringPhase < internalStep){
      alreadyReEnqueuedDuringPhase = internalStep;
      DEBUG_PRINT("WARNING: calling progress() after trailing update because work is available to process\n");
      selfContinue();
    }
    
    MERGE_PATH_DELETE_D(A,internalStep);
    
  }
  
  void processLocalLU(int ignoredParam) {
    DEBUG_PRINT("processLocalLU() called on block %d,%d\n", thisIndex.x, thisIndex.y);
    CkAssert(!done);
    
    double *incomingL = getBufferedL(internalStep);
    double *incomingU = getBufferedU(internalStep);

    CkAssert(internalStep==thisIndex.x && internalStep==thisIndex.y);

//    double mem = CmiMemoryUsage();
//    CkPrintf("CmiMemoryUsage() = %lf\n", mem);
    
    // We are the top-left-most active block
    
    // Verify that there are no outstanding buffered messages.
    CkAssert(buffersEmpty());

    // CkPrintf("[%d] chare %d,%d internalStep=%d solveLocalLU\n", CkMyPe(), thisIndex.x, thisIndex.y, internalStep);
    solveLocalLU(); // compute a local LU
      

    // If this is the very last bottom rightmost block
    if (thisIndex.x == numBlks-1 && thisIndex.y == numBlks-1) {
      DEBUG_PRINT("[%d] chare %d,%d calling mainProxy.iterationCompleted\n", CkMyPe(), thisIndex.x, thisIndex.y);
      mainProxy.iterationCompleted();
    } else {
      DEBUG_PRINT("[%d] chare %d,%d is top,left block this step, multicast U & L\n", CkMyPe(), thisIndex.x, thisIndex.y);
      multicastRecvU();	//broadcast the U downwards to the blocks in the same column
      multicastRecvL(); 	//broadcast the L rightwards to the blocks in the same row
    }
      
    MERGE_PATH_DELETE_D(A,internalStep);
    // This block is now done
    DEBUG_PRINT("chare %d,%d is now done\n",  thisIndex.x, thisIndex.y);
    done = true;
    
  }

  /// Store a pointer to the buffered L messages
  void bufferL(blkMsg *msg) {
    CmiReference(UsrToEnv(msg));
    //    DEBUG_PRINT("elem %d,%d buffering UBuffers[%d] %p\n", thisIndex.x, thisIndex.y, msg->step, msg);
    for (int i=LBuffers.size(); i<msg->step; i++) 
      LBuffers.insert(i, (blkMsg*)NULL);
    LBuffers.insert(msg->step, msg);
  }
  
  /// Store a pointer to the buffered U message
  void bufferU(blkMsg *msg) {
    CmiReference(UsrToEnv(msg));
    //   DEBUG_PRINT("elem %d,%d buffering UBuffers[%d] %p\n", thisIndex.x, thisIndex.y, msg->step, msg);
    for (int i=UBuffers.size(); i<msg->step; i++) 
      UBuffers.insert(i, (blkMsg*)NULL);
    UBuffers.insert(msg->step, msg);
  }
  
  /// Call progress on myself, possibly using priorities 
  inline void selfContinue(){
    int integerPrio;
    
    //    CkPrintf("continuing %d,%d  internalStep=%d \n", thisIndex.x,thisIndex.y, internalStep);


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

  /// Are there enough buffered messages for this step to perform the required work?
  continueWithTask canContinue(){
    double *incomingL = getBufferedL(internalStep);
    double *incomingU = getBufferedU(internalStep);
    
    if(incomingL!=NULL && incomingU!=NULL) 
      return CONTINUE_TRAIL;
    
    if (internalStep==thisIndex.x && incomingL!=NULL)
      return CONTINUE_U;
    
    if (internalStep==thisIndex.y && incomingU!=NULL)
      return CONTINUE_L;    
    
    if(internalStep==thisIndex.x && internalStep==thisIndex.y)
      return CONTINUE_LU;
    
    return NO_CONTINUE;    
  }
  
 bool buffersEmpty(){
    int validMsgInUBuffer = 0;
    int validMsgInLBuffer = 0;
    
    for (int i=0; i<UBuffers.size(); i++) 
      if(UBuffers[i] != NULL)
	validMsgInUBuffer ++;
	
    for (int i=0; i<LBuffers.size(); i++) 
      if(LBuffers[i] != NULL)
	validMsgInLBuffer ++;
    
    return validMsgInUBuffer==0 && validMsgInLBuffer==0;
  }

  void printInfo(){
    int validMsgInUBuffer = 0;
    int validMsgInLBuffer = 0;
    
    for (int i=0; i<UBuffers.size(); i++) 
      if(UBuffers[i] != NULL)
	validMsgInUBuffer ++;
	
    for (int i=0; i<LBuffers.size(); i++) 
      if(LBuffers[i] != NULL)
	validMsgInLBuffer ++;
    
    //  CkPrintf("[%d] elem %d,%d has %d buffered U messages and %d buffered L messages\n", CkMyPe(), thisIndex.x, thisIndex.y, validMsgInUBuffer, validMsgInLBuffer);
    
  }

  double *getBufferedL(int idx) {
    //    DEBUG_PRINT("getBufferedL: size=%d, idx=%d\n", LBuffers.size(), idx );
    if (LBuffers.size() <= idx) return NULL;
    if (LBuffers[idx] == NULL) return NULL;

    double *ret = LBuffers[idx]->data;
    CkAssert(ret!=NULL);
    // DEBUG_PRINT("idx=%d LBuffers[idx]->step=%d\n", idx,LBuffers[idx]->step);

    CkAssert(idx==LBuffers[idx]->step);

    return ret;
  }

  double *getBufferedU(int idx) {
    //  DEBUG_PRINT("getBufferedU: size=%d, idx=%d\n", UBuffers.size(), idx );
    if (UBuffers.size() <= idx) return NULL;
    if (UBuffers[idx] == NULL) return NULL;

    double *ret = UBuffers[idx]->data;
    CkAssert(ret!=NULL);
    CkAssert(idx==UBuffers[idx]->step);

    return ret;
  }
  
  void removeBufferedL(int idx) {
    //  DEBUG_PRINT("elem %d,%d deleting LBuffers[%d] %p\n", thisIndex.x, thisIndex.y, idx, LBuffers[idx]);
    CkAssert(idx < LBuffers.size());
    CmiFree(UsrToEnv(LBuffers[idx]));
    LBuffers[idx]=NULL;
  }
  
  void removeBufferedU(int idx) {
    //  DEBUG_PRINT("elem %d,%d deleting UBuffers[%d] %p\n", thisIndex.x, thisIndex.y, idx, UBuffers[idx]);
    CkAssert(idx < UBuffers.size());
    CmiFree(UsrToEnv(UBuffers[idx]));
    UBuffers[idx]=NULL;
  }


  void deallocateBuffers() {
    UBuffers.free();
    LBuffers.free();
    MERGE_PATH_DELETE_ALL_D(A);
  }

private:
  //internal functions for creating messages to encapsulate the priority
  inline blkMsg* createABlkMsg() {
    blkMsg *msg;
    
    if(doPrioritize) {
      msg = new(BLKSIZE*BLKSIZE, 8*sizeof(int)) blkMsg(BLKSIZE);
      DEBUG_PRINT("setting priority to internalStep=%d\n", internalStep);
      *((int*)CkPriorityPtr(msg)) = (int)internalStep;
      CkSetQueueing(msg, CK_QUEUEING_IFIFO);
    } else {
      msg = new(BLKSIZE*BLKSIZE)blkMsg(BLKSIZE);
    }
    
    return msg;
  }

  
};

#include "lu.def.h"
