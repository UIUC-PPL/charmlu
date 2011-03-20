#ifndef LU_PARAMS
#define LU_PARAMS

#include "lu.decl.h"

// Internal configuration and parameters for LUSolver
struct LUParams {
  // Trace handles for identifying regions of code
  int traceTrailingUpdate, traceComputeU, traceComputeL, traceSolveLocalLU;
  CkGroupID mcastMgrGID;

  void pup(PUP::er &p) {
    p | traceTrailingUpdate;
    p | traceComputeU;
    p | traceComputeL;
    p | traceSolveLocalLU;
    p | mcastMgrGID;
  }
};

#endif
