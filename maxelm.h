#ifndef CHARMLU_MAXELMH
#define CHARMLU_MAXELMH

#include <charm++.h>

struct MaxElm {
  double val;
  int loc;
  MaxElm(): val(0.0), loc(-1) { }
  MaxElm(double _val, int _loc): val(_val), loc(_loc) { }
};
PUPbytes(MaxElm)

#endif
