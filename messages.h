#ifndef LU_MESSAGES_H
#define LU_MESSAGES_H

#include "lu.decl.h"
#include <algorithm>

/// A block of the matrix, with ornamentation for setup-free multicast
struct blkMsg: public CkMcastBaseMsg, CMessage_blkMsg {
  /// Payload
  double *data;
  CkIndex2D src;

  /// Destinations
  int *pes;
  int npes_sender, npes_receiver, offset;
  /// Whether the message has gone to the first half of its destination list
  bool firstHalfSent;

  blkMsg(CkIndex2D index_)
    : src(index_), npes_sender(0), npes_receiver(0), offset(0),
      firstHalfSent(false) {
    CkSetRefNum(this, std::min(src.x, src.y));
    CkSetQueueing(this, CK_QUEUEING_IFIFO);
  }
};

struct BlockReadyMsg : public CkMcastBaseMsg, CMessage_BlockReadyMsg {
  BlockReadyMsg(CkIndex2D idx) : src(idx)
  {
    CkSetRefNum(this, std::min(src.x, src.y));
  }
  CkIndex2D src;
};

class rednSetupMsg: public CkMcastBaseMsg, public CMessage_rednSetupMsg
{
    public:
        CkGroupID rednMgrGID;
        rednSetupMsg(CkGroupID _gid): rednMgrGID(_gid) {}
};

class BVecMsg : public CMessage_BVecMsg, public CkMcastBaseMsg {
public:
  BVecMsg(const int numElements, double *bvec, bool _forward) {
    memcpy(data, bvec, numElements * sizeof(double));
    forward = _forward;
  }
  /// The post-diagonal portion of the active row
  double *data;
  bool forward;
};

class UMsg : public CMessage_UMsg, public CkMcastBaseMsg {
    public:
        UMsg(const int numElements, double *useg) {
            memcpy(data, useg, numElements * sizeof(double));
        }
        /// The post-diagonal portion of the active row
        double *data;
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

#endif
