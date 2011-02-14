
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
