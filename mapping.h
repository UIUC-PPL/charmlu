#ifndef LU_MAPPING_H
#define LU_MAPPING_H

#include "lu.decl.h"
#include <string>
#include <set>

// Implement a mapping that tiles a 2D processor tile in the 2D chare array
struct PE2DTilingMap : public CBase_PE2DTilingMap {
  const int peRows, peCols, peRotate, peStride, numBlks;

  PE2DTilingMap(int _peRows, int _peCols, int _peRotate, int _peStride, int _numBlks)
    : peRows(_peRows), peCols(_peCols), peRotate(_peRotate), peStride(_peStride), numBlks(_numBlks) { }

  int map(const int coor[2]) {
    int tileYIndex = coor[1]  / peCols;
    int XwithinPEtile = (coor[0] + tileYIndex * peRotate) % peRows;
    int YwithinPEtile = coor[1] % (peCols / peStride);
    int subtileY = (coor[1] % peCols) / (peCols / peStride);
    int peNum = XwithinPEtile * peStride + YwithinPEtile * peStride * peRows + subtileY;
    return peNum;
  }

  int procNum(int arrayHdl, const CkArrayIndex &idx) { return map((int*)idx.data()); }

  int pesInPanel(CkIndex2D index) {
    if (index.x < index.y) return peRows;
    else {
      std::set<int> pesInRow;
      int coor[2];

      coor[0] = index.x;
      for (int y = index.y+1; y < numBlks; ++y) {
	coor[1] = y;
	pesInRow.insert(map(coor));
      }

      // Leave out the PE that the block lives on - it will never send to itself
      coor[0] = index.x; coor[1] = index.y;
      pesInRow.erase(map(coor));

      return pesInRow.size();
    }
  }
};

/// Mapping for BlockScheduler - created as an array so group section multicast
/// can be used
struct OnePerPE : public CBase_OnePerPE {
  int registerArray(CkArrayIndexMax& numElements,CkArrayID aid) { return 0; }
  int procNum(int arrayHdl, const CkArrayIndex &elt) { return *(int*)elt.data(); }
};

#endif /* LU_MAPPING_H */
