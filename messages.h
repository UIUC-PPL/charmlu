#ifndef LU_MESSAGES_H
#define LU_MESSAGES_H

#include "luMessages.decl.h"
#include <ckmulticast.h>
#include <algorithm>

/// Utility to retain a reference to a reference-counted message
static inline void takeRef(void *m) { CmiReference(UsrToEnv(m)); }
/// Utility to drop a reference to a reference-counted message
static inline void dropRef(void *m) { CmiFree(UsrToEnv(m)); }

/// A block of the matrix, with ornamentation for setup-free multicast
struct blkMsg: public CkMcastBaseMsg, CMessage_blkMsg {
  double *data;

  blkMsg(double *data_, int step, int blkSize) {
    CkSetRefNum(this, step);
    memcpy(data, data_, sizeof(double) * blkSize * blkSize);
    //CkSetQueueing(this, CK_QUEUEING_IFIFO);
  }
};

struct rednSetupMsg: public CkMcastBaseMsg, public CMessage_rednSetupMsg {
  CkGroupID rednMgrGID;

  rednSetupMsg(CkGroupID _gid): rednMgrGID(_gid) {}
};

struct BVecMsg : public CMessage_BVecMsg, public CkMcastBaseMsg {
  /// The post-diagonal portion of the active row
  double *data;
  bool forward;

  BVecMsg(const int numElements, double *bvec, bool _forward) {
    memcpy(data, bvec, numElements * sizeof(double));
    forward = _forward;
  }
};

#endif
