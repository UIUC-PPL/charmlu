#include "lu.decl.h"
#include "lu.h"

void BlockScheduler::registerBlock(CkIndex2D index) {

}

void BlockScheduler::pivotsDone(CkIndex2D index) {

}

void BlockScheduler::dataReady(CkIndex2D index, BlockReadyMsg *m) {

}

inline bool operator==(const CkIndex2D &l, const CkIndex2D &r)
{ return l.x == r.x && l.y == r.y; }
inline bool operator<(const CkIndex2D &l, const CkIndex2D &r)
{ return l.x < r.x || (l.x == r.x && l.y < r.y); }

