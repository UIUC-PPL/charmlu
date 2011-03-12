#ifndef LU_MANAGER_H
#define LU_MANAGER_H

#include "lu.decl.h"
#include "messages.h"

// Prefix notation:
// b: an index into an LUBlk chare array

enum PrioType {
  SEND_PIVOT_DATA, PROCESS_TRAILING_UPDATE, RECVL, DIAG_SEND_PIVOT,
  BELOW_SEND_USEG, DIAG_SEND_USEG, MULT_RECV_U, MULT_RECV_L,
  PIVOT_RIGHT_SEC, PIVOT_LEFT_SEC, PIVOT_CRITICAL, PIVOT_NOT_CRITICAL,
  GET_BLOCK, SEND_BLOCKS, DIAG_MULT_RECV_L
};

class LUMgr : public CBase_LUMgr
{
protected:
  int BLKSIZE;
  int matSize;
  int numBlks;
  LUMgr(int BLKSIZE_, int matSize_) : BLKSIZE(BLKSIZE_), matSize(matSize_) {
    numBlks = matSize / BLKSIZE;
  }

public:
  /// Return a message with space for a block's data and the requested
  /// priority bits. Subclasses may actually set the priorities
  /// appropriately.
  virtual blkMsg *createBlockMessage(int brow, int bcol, int bactive,
                                     int priobits = 0)
  {
    return new (BLKSIZE*BLKSIZE, priobits) blkMsg;
  }
  virtual CkEntryOptions& setPrio(PrioType type, CkEntryOptions& opts,
                                  int y = 0 , int x = 0) = 0;
  virtual void setPrio(CkMessage *msg, PrioType type, int refnum = -1, int y = 0, int x = 0) = 0;
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

  CkEntryOptions& setPrio(PrioType type, CkEntryOptions& opts, int y = 0, int x = 0) {
    int prio = 0;
    switch (type) {
    case MULT_RECV_L:
      prio = x * BLKSIZE;
      break;
    case MULT_RECV_U:
      prio = (y - 3) * BLKSIZE;
      break;
    case SEND_PIVOT_DATA: case DIAG_SEND_PIVOT:
      prio = -1;
      break;
    case RECVL:
      prio = 1;
      break;
    case GET_BLOCK:
      prio = -3;
      break;
    case SEND_BLOCKS:
      prio = -2;
      break;
    case PROCESS_TRAILING_UPDATE:
      prio = y * BLKSIZE;
      break;
    }
    opts.setPriority(prio);
    return opts;
  }

  void setPrio(CkMessage *msg, PrioType type, int refnum = -1, int y = 0, int x = 0) {
    int prio = 0;
    switch (type) {
    case DIAG_SEND_USEG: case DIAG_MULT_RECV_L: case BELOW_SEND_USEG:
      prio = -1;
      break;
    case PIVOT_RIGHT_SEC:
      prio = y * BLKSIZE;
      break;
    case PIVOT_LEFT_SEC:
      prio = BLKSIZE * numBlks + 1;
      break;
    case PIVOT_CRITICAL:
      prio = (y - 1) * BLKSIZE;
      break;
    case PIVOT_NOT_CRITICAL:
      prio = numBlks * BLKSIZE;
      break;
    }
    *(int*)CkPriorityPtr(msg) = prio;
    CkSetQueueing(msg, CK_QUEUEING_IFIFO);
    // Better way than this?
    if (refnum != -1)
      CkSetRefNum(msg, refnum);
  }
};

#endif
