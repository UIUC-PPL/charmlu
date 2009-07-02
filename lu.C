#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <pthread.h>
 
#include <malloc.h>


#include <comlib.h>
#include <controlPoints.h> // must come before user decl.h if they are using the pathInformationMsg
#include "lu.decl.h"

#if USE_CBLAS_H
extern "C" {
#include <cblas.h>
#include <clapack.h>
}
#elif USE_MKL_CBLAS_H
#include "mkl_cblas.h"
#include "mkl_lapack.h"
#elif USE_ACCELERATE_BLAS
#include <Accelerate/Accelerate.h>
#else
#error "No BLAS Header files included!"
#endif


#define numIterations 1


/* readonly: */
CProxy_Main mainProxy;
CProxy_LUBlk luArrProxy;
ComlibInstanceHandle cinst0; 
ComlibInstanceHandle cinst1; 
ComlibInstanceHandle cinst2;  
int gMatSize;
int numBlks;
int BLKSIZE;
int whichMapping;
int traceTrailingUpdate;
int traceComputeU;
int traceComputeL;
int traceSolveLocalLU;
int doPrioritize;

 
//#define DEBUG_PRINT(...) CkPrintf(__VA_ARGS__)
#define DEBUG_PRINT(...) 


/// Called periodically for debugging
static void periodicDebug(void* ptr, double currWallTime){
  //  CcdCallFnAfterOnPE((CcdVoidFn)periodicDebug, (void*)NULL, controlPointSamplePeriod, CkMyPe());
  luArrProxy.printInfo();
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

  void getNRndInts(int num, int *d) {
    for (int i=0; i<num; i++)
      d[i] = nextRndInt();
  }
};


class blkMsg: public CMessage_blkMsg {
public:
  double *data;
  int step;

  void setMsgData(double *d, int s) {
    memcpy(data, d, sizeof(double)*BLKSIZE*BLKSIZE);
    step = s;
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

  Main(CkArgMsg* m) {
    iteration = 0;

    if (m->argc<2) {
      CkPrintf("Usage: %s <matrix size>\n", m->argv[0]);
      CkExit();
    }
    
    mainProxy = thisProxy;
    
      
    DEBUG_PRINT("Registering user events\n");
    traceTrailingUpdate = traceRegisterUserEvent("Trailing Update");
    traceComputeU = traceRegisterUserEvent("Compute U");
    traceComputeL = traceRegisterUserEvent("Compute L");
    traceSolveLocalLU = traceRegisterUserEvent("Solve local LU");

    traceRegisterUserEvent("Local Multicast Deliveries", 10000);    
    traceRegisterUserEvent("Remote Multicast Forwarding", 10001);

    BLKSIZE = 1 << staticPoint("Block Size", 8,8);
    whichMapping = 0;
    doPrioritize = 0;

    gMatSize = atoi(m->argv[1]);
 
    if (gMatSize%BLKSIZE!=0) 
      CkAbort("The matrix size should be a multiple of block size!\n");
  
    numBlks = gMatSize/BLKSIZE;
  
    CkPrintf("Running LU on %d processors (%d nodes) on matrix %dX%d with block size %d\n",
	     CkNumPes(), CmiNumNodes(), gMatSize, gMatSize, BLKSIZE);
    
 //    Strategy * strategy0 = new DirectMulticastStrategy();
//     cinst0 = ComlibRegister(strategy0);
    
//     Strategy * strategy1 = new RingMulticastStrategy();
//     cinst1 = ComlibRegister(strategy1);
    
    Strategy * strategy2 = new OneTimeNodeTreeMulticastStrategy();     
    cinst2 = ComlibRegister(strategy2); 
    
    if(whichMapping==0){
      CProxy_BlockCyclicMap myMap = CProxy_BlockCyclicMap::ckNew();
      CkArrayOptions opts(numBlks, numBlks);
      opts.setMap(myMap);
      luArrProxy = CProxy_LUBlk::ckNew(opts);
    } else {
      luArrProxy = CProxy_LUBlk::ckNew(numBlks, numBlks);
    }
    
    CkCallback *cb = new CkCallback(CkIndex_Main::arrayIsCreated(NULL), thisProxy);
    CkStartQD(*cb); // required currently for use with new Comlib
    
    CcdCallFnAfterOnPE((CcdVoidFn)periodicDebug, (void*)NULL, 10000, CkMyPe());
    
  }

  
  void arrayIsCreated(CkReductionMsg * m) {
    startTime = CmiWallTimer();
    luArrProxy(0,0).begin();
  }
  

  
  void iterationCompleted() {
    double endTime = CmiWallTimer();
    double duration = endTime-startTime;
    registerControlPointTiming(duration); 
    
    
    CkPrintf("Iteration %d time: %fs\n", iteration, duration);
    if(iteration == numIterations-1){ 
      terminateProg();
      return;
    }
    
    iteration++;
    
    //    gotoNextPhase();
    traceUserSuppliedNote("*** New Iteration");
    
    staticPoint("Block Size", 8,8); // call this so it gets recorded for this phase
    int whichMulticastStrategy = controlPoint("which multicast strategy", 3,3);
    
    CkCallback *cb = new CkCallback(CkIndex_Main::arrayIsCreated(NULL), thisProxy); 
    luArrProxy.ckSetReductionClient(cb);
    luArrProxy.init(whichMulticastStrategy);
    
  }
  
  
  void terminateProg() {
    double endTime = CmiWallTimer();
    double duration = endTime-startTime;

    CkPrintf("Main execution time: %fs\n", duration);
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

    std::cout << "If run on order.cs.uiuc.edu, I think you are getting  \t" << fractionOfPeakOnOrderPercent << "% of peak" << std::endl;


    registerControlPointTiming(duration);

    CkCallback cb(CkIndex_Main::done(NULL),thisProxy); 
    traceCriticalPathBack(cb);

  }


  void done(pathInformationMsg *m){
    // CkPrintf("Main::done() After critical path has been determined\n");
    //    m->printme();
    CkExit();
  }

};


class LUBlk: public CBase_LUBlk {

public:
  
  double *LU;
  int whichMulticastStrategy;
  bool done;

  int alreadyReEnqueuedDuringPhase;
  
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

    CkAssert(BLKSIZE>0); // If this fails, readonly variables aren't
			 // propagated soon enough. I'm assuming they
			 // are safe to use here.

    //    LU = new double[BLKSIZE*BLKSIZE];
    LU = (double*)memalign(128, BLKSIZE*BLKSIZE*sizeof(double) );
    //   CkPrintf("LU mod 128 = %lu\n", ((unsigned long)LU) % 128);


    internalStep = 0;  
     
    traceUserSuppliedData(-1);  
    traceMemoryUsage();  
     
    MatGen rnd(thisIndex.x * numBlks + thisIndex.y);      
    for (int i=0; i<BLKSIZE*BLKSIZE; i++) {  
      LU[i] = rnd.toRndDouble(rnd.nextRndInt());  
    }  

  }


  void init(int _whichMulticastStrategy){
    
    done = false;
    alreadyReEnqueuedDuringPhase = -1;

    whichMulticastStrategy = _whichMulticastStrategy;
    
    internalStep = 0; 
    
    traceUserSuppliedData(-1); 
    traceMemoryUsage(); 
    
    MatGen rnd(thisIndex.x * numBlks + thisIndex.y);     
    for (int i=0; i<BLKSIZE*BLKSIZE; i++) { 
      LU[i] = rnd.toRndDouble(rnd.nextRndInt()); 
    } 

    contribute();
  }
  

  ~LUBlk() {

    delete [] LU;
    LU = NULL;

  }


  LUBlk(CkMigrateMessage* m) {}

  //added for migration
  void pup(PUP::er &p) {  }


 
  // Called on element 0,0
  void begin() {
    CkAssert(thisIndex.x == 0 && thisIndex.y == 0);
    traceUserSuppliedData(internalStep);
    traceMemoryUsage();
    CkAssert(!done);
    thisProxy(0,0).progress();
  }


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
#warning "Using lapack for solveLocalLU"
    //There's an output for permuation array which is not
    //used in the current non-lapack version. This permuation
    //array should also be broadcasted to those elements that
    //only needs ComputeL and ComputeU (???)
    int *ipiv = new int[BLKSIZE];
    
#if 1
    clapack_dgetrf(CblasRowMajor, BLKSIZE, BLKSIZE, LU, BLKSIZE, ipiv);
#else
    int info;
    // This one doesn't quite do what we want... it does pivoting
    dgetrf(&BLKSIZE, &BLKSIZE, LU, &BLKSIZE, ipiv, &info);
#endif
    
    
    //  CkPrintf("Permutation info: ");
//     for (int i=0; i<BLKSIZE; i++) CkPrintf("%d, ", ipiv[i]);
//     CkPrintf("\n");
    delete [] ipiv;
// #else
//     //first compute local LU
//     for (int iter=0; iter<BLKSIZE; iter++) {
//       double pivot = LU[iter*BLKSIZE+iter];
//       //compute L of [iter]th column
//       //L[iter*BLKSIZE+iter]=1.0;
//       for (int row=iter+1; row<BLKSIZE; row++) {
// 	LU[row*BLKSIZE+iter] = LU[row*BLKSIZE+iter]/pivot;
// 	//U[row*BLKSIZE+iter] = 0.0;
//       }

//       //The [iter]th row of U need not to be changed

//       //update the remaining matrix
//       for (int row=iter+1; row<BLKSIZE; row++)
// 	for (int col=iter+1; col<BLKSIZE; col++)
// 	  LU[row*BLKSIZE+col] -= LU[row*BLKSIZE+iter]*LU[iter*BLKSIZE+col];
//     }
// #endif

    traceUserBracketEvent(traceSolveLocalLU, luStart, CmiWallTimer());
  }

  // really is compute W
  void computeU(blkMsg *givenLMsg) {
    traceUserSuppliedData(internalStep);
    traceMemoryUsage();
    double *givenL = givenLMsg->data;

    double uStart = CmiWallTimer();

    DEBUG_PRINT("elem[%d,%d]::computeU called at step %d\n", thisIndex.x, thisIndex.y, internalStep);

    //processing row by row (forward substitution)
    //the 1st row of U is not changed

    //solve following rows based on previously solved rows
    //row indicates the row of U that is just solved

    cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, BLKSIZE, BLKSIZE, 1.0, givenL, BLKSIZE, LU, BLKSIZE);
    

    traceUserBracketEvent(traceComputeU, uStart, CmiWallTimer());
  }

  void computeL(blkMsg *givenUMsg) {
    traceUserSuppliedData(internalStep);
    traceMemoryUsage();
    double *givenU = givenUMsg->data;


    //traceUserEvent(traceComputeL);
    double lStart = CmiWallTimer();

    DEBUG_PRINT("elem[%d,%d]::computeL called at step %d\n", thisIndex.x, thisIndex.y, internalStep);


    cblas_dtrsm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, BLKSIZE, BLKSIZE, 1.0, givenU, BLKSIZE, LU, BLKSIZE);
    
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

    cblas_dgemm( CblasRowMajor,
		 CblasNoTrans, CblasNoTrans,
		 BLKSIZE, BLKSIZE, BLKSIZE,
		 -1.0, incomingL,
		 BLKSIZE, incomingU, BLKSIZE,
		 1.0, LU, BLKSIZE);
    

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
      thisProxy(thisIndex.x,thisIndex.y).progress();
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
      thisProxy(thisIndex.x,thisIndex.y).progress();
    }
  }

  //broadcast the U downwards to the blocks in the same column
  inline void multicastRecvU() {
    traceUserSuppliedData(internalStep);
    traceMemoryUsage();
    
    DEBUG_PRINT("[PE %d] elem %d,%d Multicast to part of column %d from step %d\n", CkMyPe(), thisIndex.x, thisIndex.y, thisIndex.y, internalStep);
    
    CProxySection_LUBlk oneCol = CProxySection_LUBlk::ckNew(thisArrayID, thisIndex.x+1, numBlks-1, 1, thisIndex.y, thisIndex.y, 1);
    

  //   switch(whichMulticastStrategy){
//     case 0:
//       // no delegation
//       break;
//     case 1:
//       CkAssert(cinst0);
//       ComlibAssociateProxy(cinst0, oneCol);
//       break;
//     case 2: 
//       CkAssert(cinst1); 
//       ComlibAssociateProxy(cinst1, oneCol);        
//       break;
//     case 3:
      CkAssert(cinst2);
      ComlibAssociateProxy(cinst2, oneCol);
    //   break;
//     }
    
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
    
//     switch(whichMulticastStrategy){ 
//     case 0: 
//       // no delegation 
//       break;
//     case 1:
//       CkAssert(cinst0);
//       ComlibAssociateProxy(cinst0, oneRow); 
//       break; 
//     case 2:  
//       CkAssert(cinst1);  
//       ComlibAssociateProxy(cinst1, oneRow);         
//       break;
//     case 3:  
      CkAssert(cinst2);  
      ComlibAssociateProxy(cinst2, oneRow);         
//       break;
//     }
    
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
  
 

  void progress() {
    DEBUG_PRINT("progress() called on block %d,%d\n", thisIndex.x, thisIndex.y);
    CkAssert(!done);
    
    double *incomingL = getBufferedL(internalStep);
    double *incomingU = getBufferedU(internalStep);

    if( internalStep==thisIndex.x && internalStep==thisIndex.y ) {
      // If we are the top-left-most active block
      	
      // Verify that there are no outstanding buffered messages.
      CkAssert(buffersEmpty());

      //      CkPrintf("[%d] chare %d,%d internalStep=%d solveLocalLU\n", CkMyPe(), thisIndex.x, thisIndex.y, internalStep);
      solveLocalLU(); // compute a local LU
      

      // If this is the very last bottom rightmost block
      if (thisIndex.x == numBlks-1 && thisIndex.y == numBlks-1) {
	DEBUG_PRINT("[%d] chare %d,%d calling mainProxy.iterationCompleted -------------------------------------------------\n", CkMyPe(), thisIndex.x, thisIndex.y);
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
    else if (internalStep==thisIndex.x && incomingL!=NULL) {
      // Else, if we are in the top row of active blocks, and we
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
    else if (internalStep==thisIndex.y && incomingU!=NULL) {
      // Else, if we are in the left row of active blocks
      
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
    else if(incomingL!=NULL && incomingU!=NULL) {
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
	DEBUG_PRINT("WARNING:                                  calling progress() after trailing update because work is available to process\n");
	thisProxy(thisIndex.x,thisIndex.y).progress();
      }
      
      MERGE_PATH_DELETE_D(A,internalStep);
    }
    
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
  

  /// Are there enough buffered messages for this step to perform the required work 
  bool canContinue(){
    double *incomingL = getBufferedL(internalStep);
    double *incomingU = getBufferedU(internalStep);
    
    if(incomingL!=NULL && incomingU!=NULL) 
      return true;
    
    if (internalStep==thisIndex.x && incomingL!=NULL)
      return true;
    
    if (internalStep==thisIndex.y && incomingU!=NULL)
      return true;    
    
    if(internalStep==thisIndex.x && internalStep==thisIndex.y)
      return true;
    
    return false;    
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
  //internal functions for creating messages to encapsulate the priority function of functions
  inline blkMsg* createABlkMsg() {
    blkMsg *msg;
    
    if(doPrioritize) {
      msg = new(BLKSIZE*BLKSIZE, 8*sizeof(int)) blkMsg;
      DEBUG_PRINT("setting priority to internalStep=%d\n", internalStep);
      *((int*)CkPriorityPtr(msg)) = (int)internalStep;
      CkSetQueueing(msg, CK_QUEUEING_IFIFO);
    } else {
      msg = new(BLKSIZE*BLKSIZE)blkMsg;
    }
    
    return msg;
  }
  

};


#include "lu.def.h"
