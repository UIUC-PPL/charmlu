
// Prefix notation:
// b: an index into an LUBlk chare array

class LUMgr : public CBase_LUMgr
{
protected:
  int BLKSIZE;
  int matSize;
  LUMgr(int BLKSIZE_, int matSize_) : BLKSIZE(BLKSIZE_), matSize(matSize_) {}

public:
  /// Return a message with space for a block's data and the requested
  /// priority bits. Subclasses may actually set the priorities
  /// appropriately.
  virtual blkMsg *createBlockMessage(int brow, int bcol, int bactive,
                                     int priobits = 0)
  {
    return new (BLKSIZE*BLKSIZE, priobits) blkMsg;
  }
};

struct PrioLU : public LUMgr
{
  PrioLU(int BLKSIZE, int matSize) : LUMgr(BLKSIZE, matSize) {}

  blkMsg *createBlockMessage(int brow, int bcol, int bactive, int priobits)
  {
    blkMsg *msg = LUMgr::createBlockMessage(brow, bcol, bactive, priobits);
    CkSetQueueing(msg, CK_QUEUEING_IFIFO);
    if (0 != priobits)
      *((int*)CkPriorityPtr(msg)) = bactive;
    return msg;
  }
};


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
