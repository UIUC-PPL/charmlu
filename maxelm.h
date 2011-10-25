#ifndef CHARMLU_MAXELMH
#define CHARMLU_MAXELMH

#include <charm++.h>

/// A pair representing the value and the location of the candidate pivot element
struct MaxElm {
  double val;
  int loc;
  MaxElm(double _val=0.0, int _loc=-1): val(_val), loc(_loc) { }
};
PUPbytes(MaxElm)

/// The ID of a custom reduction function for MaxElm types
extern CkReduction::reducerType MaxElmReducer;

#endif
