#ifndef LU_MESSAGES_H
#define LU_MESSAGES_H

#include "lu.decl.h"
#include <algorithm>

struct blkMsg: public CkMcastBaseMsg, CMessage_blkMsg {
  double *data;
  CkArrayIndex1D *pes;
  int npes, offset;
  CkIndex2D src;

  void setMsgData(double *data_, int step, int BLKSIZE, CkIndex2D index) {
    memcpy(data, data_, sizeof(double)*BLKSIZE*BLKSIZE);
    src = index;
    CkSetRefNum(this, step);
  }
};

void propagateBlkMsg(blkMsg *msg, CProxy_BlockScheduler bs);

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
