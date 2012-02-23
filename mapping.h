#ifndef LU_MAPPING_H
#define LU_MAPPING_H

#include "lu.decl.h"
#include <string>

struct LUMap : public CkArrayMap {
  LUMap() {}
  virtual std::string desc() { return ""; }
  virtual int pesInPanel(CkIndex2D index) { return CkNumPes(); }
  int procNum(int arrayHdl, const CkArrayIndex &idx) { return 0; }
};

class BlockCyclicMap : public LUMap {
  int r, num_blocks;
public:
  BlockCyclicMap(int r_, int num_blocks_) : r(r_), num_blocks(num_blocks_) {}
  int procNum(int arrayHdl, const CkArrayIndex &idx) {
    int *coor = (int *)idx.data();
    int m = coor[1] * num_blocks + coor[0];
    return (m % (r * CkNumPes())) / r;
  }
};

#include <set>

// Implement a mapping that tiles a 2D processor tile in the 2D chare array
class PE2DTilingMap: public LUMap {
public:
  PE2DTilingMap(int _peRows, int _peCols, int _peRotate, int _peStride, int _numBlks)
    : peRows(_peRows), peCols(_peCols), peRotate(_peRotate), peStride(_peStride),
      numBlks(_numBlks) {
    CkAssert(peRows > 0 && peCols > 0);
  }

  int map(const int coor[2]) {


if (coor[0] == 11 && coor[1] == 11) return 29;
else if (coor[0] == 11 && coor[1] == 7) return 62;
else if (coor[0] == 11 && coor[1] == 17) return 12;
else if (coor[0] == 11 && coor[1] == 2) return 38;
else if (coor[0] == 11 && coor[1] == 1) return 56;
else if (coor[0] == 11 && coor[1] == 18) return 40;
else if (coor[0] == 11 && coor[1] == 0) return 46;
else if (coor[0] == 11 && coor[1] == 13) return 41;
else if (coor[0] == 11 && coor[1] == 16) return 39;
else if (coor[0] == 11 && coor[1] == 6) return 36;
else if (coor[0] == 11 && coor[1] == 3) return 54;
else if (coor[0] == 11 && coor[1] == 9) return 40;
else if (coor[0] == 11 && coor[1] == 12) return 8;
else if (coor[0] == 11 && coor[1] == 14) return 7;
else if (coor[0] == 11 && coor[1] == 15) return 50;
else if (coor[0] == 11 && coor[1] == 8) return 60;
else if (coor[0] == 11 && coor[1] == 4) return 59;
else if (coor[0] == 11 && coor[1] == 19) return 54;
else if (coor[0] == 11 && coor[1] == 10) return 37;
else if (coor[0] == 11 && coor[1] == 5) return 53;
else if (coor[0] == 7 && coor[1] == 11) return 53;
else if (coor[0] == 7 && coor[1] == 7) return 62;
else if (coor[0] == 7 && coor[1] == 2) return 31;
else if (coor[0] == 7 && coor[1] == 17) return 19;
else if (coor[0] == 7 && coor[1] == 1) return 22;
else if (coor[0] == 7 && coor[1] == 18) return 18;
else if (coor[0] == 7 && coor[1] == 0) return 16;
else if (coor[0] == 7 && coor[1] == 16) return 61;
else if (coor[0] == 7 && coor[1] == 13) return 35;
else if (coor[0] == 7 && coor[1] == 6) return 50;
else if (coor[0] == 7 && coor[1] == 3) return 41;
else if (coor[0] == 7 && coor[1] == 9) return 50;
else if (coor[0] == 7 && coor[1] == 12) return 42;
else if (coor[0] == 7 && coor[1] == 14) return 23;
else if (coor[0] == 7 && coor[1] == 15) return 43;
else if (coor[0] == 7 && coor[1] == 8) return 27;
else if (coor[0] == 7 && coor[1] == 4) return 52;
else if (coor[0] == 7 && coor[1] == 10) return 55;
else if (coor[0] == 7 && coor[1] == 19) return 23;
else if (coor[0] == 7 && coor[1] == 5) return 27;
else if (coor[0] == 2 && coor[1] == 11) return 35;
else if (coor[0] == 2 && coor[1] == 7) return 30;
else if (coor[0] == 2 && coor[1] == 17) return 46;
else if (coor[0] == 2 && coor[1] == 2) return 4;
else if (coor[0] == 2 && coor[1] == 1) return 3;
else if (coor[0] == 2 && coor[1] == 18) return 44;
else if (coor[0] == 2 && coor[1] == 0) return 3;
else if (coor[0] == 2 && coor[1] == 13) return 54;
else if (coor[0] == 2 && coor[1] == 16) return 42;
else if (coor[0] == 2 && coor[1] == 6) return 23;
else if (coor[0] == 2 && coor[1] == 3) return 8;
else if (coor[0] == 2 && coor[1] == 9) return 47;
else if (coor[0] == 2 && coor[1] == 12) return 51;
else if (coor[0] == 2 && coor[1] == 14) return 30;
else if (coor[0] == 2 && coor[1] == 15) return 54;
else if (coor[0] == 2 && coor[1] == 8) return 38;
else if (coor[0] == 2 && coor[1] == 4) return 12;
else if (coor[0] == 2 && coor[1] == 10) return 57;
else if (coor[0] == 2 && coor[1] == 19) return 38;
else if (coor[0] == 2 && coor[1] == 5) return 17;
else if (coor[0] == 17 && coor[1] == 11) return 13;
else if (coor[0] == 17 && coor[1] == 7) return 20;
else if (coor[0] == 17 && coor[1] == 17) return 17;
else if (coor[0] == 17 && coor[1] == 2) return 55;
else if (coor[0] == 17 && coor[1] == 1) return 49;
else if (coor[0] == 17 && coor[1] == 18) return 38;
else if (coor[0] == 17 && coor[1] == 0) return 63;
else if (coor[0] == 17 && coor[1] == 13) return 59;
else if (coor[0] == 17 && coor[1] == 16) return 23;
else if (coor[0] == 17 && coor[1] == 6) return 43;
else if (coor[0] == 17 && coor[1] == 3) return 28;
else if (coor[0] == 17 && coor[1] == 9) return 22;
else if (coor[0] == 17 && coor[1] == 12) return 51;
else if (coor[0] == 17 && coor[1] == 15) return 3;
else if (coor[0] == 17 && coor[1] == 14) return 0;
else if (coor[0] == 17 && coor[1] == 8) return 50;
else if (coor[0] == 17 && coor[1] == 4) return 45;
else if (coor[0] == 17 && coor[1] == 19) return 58;
else if (coor[0] == 17 && coor[1] == 10) return 54;
else if (coor[0] == 17 && coor[1] == 5) return 39;
else if (coor[0] == 1 && coor[1] == 11) return 55;
else if (coor[0] == 1 && coor[1] == 7) return 21;
else if (coor[0] == 1 && coor[1] == 17) return 45;
else if (coor[0] == 1 && coor[1] == 2) return 0;
else if (coor[0] == 1 && coor[1] == 1) return 1;
else if (coor[0] == 1 && coor[1] == 18) return 58;
else if (coor[0] == 1 && coor[1] == 0) return 1;
else if (coor[0] == 1 && coor[1] == 16) return 47;
else if (coor[0] == 1 && coor[1] == 13) return 37;
else if (coor[0] == 1 && coor[1] == 6) return 15;
else if (coor[0] == 1 && coor[1] == 3) return 2;
else if (coor[0] == 1 && coor[1] == 9) return 36;
else if (coor[0] == 1 && coor[1] == 12) return 33;
else if (coor[0] == 1 && coor[1] == 14) return 50;
else if (coor[0] == 1 && coor[1] == 15) return 54;
else if (coor[0] == 1 && coor[1] == 8) return 28;
else if (coor[0] == 1 && coor[1] == 4) return 6;
else if (coor[0] == 1 && coor[1] == 19) return 53;
else if (coor[0] == 1 && coor[1] == 10) return 45;
else if (coor[0] == 1 && coor[1] == 5) return 10;
else if (coor[0] == 18 && coor[1] == 11) return 44;
else if (coor[0] == 18 && coor[1] == 7) return 58;
else if (coor[0] == 18 && coor[1] == 17) return 41;
else if (coor[0] == 18 && coor[1] == 2) return 59;
else if (coor[0] == 18 && coor[1] == 1) return 36;
else if (coor[0] == 18 && coor[1] == 18) return 4;
else if (coor[0] == 18 && coor[1] == 0) return 50;
else if (coor[0] == 18 && coor[1] == 13) return 49;
else if (coor[0] == 18 && coor[1] == 16) return 46;
else if (coor[0] == 18 && coor[1] == 6) return 25;
else if (coor[0] == 18 && coor[1] == 3) return 28;
else if (coor[0] == 18 && coor[1] == 9) return 17;
else if (coor[0] == 18 && coor[1] == 12) return 1;
else if (coor[0] == 18 && coor[1] == 14) return 24;
else if (coor[0] == 18 && coor[1] == 15) return 49;
else if (coor[0] == 18 && coor[1] == 8) return 39;
else if (coor[0] == 18 && coor[1] == 4) return 44;
else if (coor[0] == 18 && coor[1] == 19) return 0;
else if (coor[0] == 18 && coor[1] == 10) return 8;
else if (coor[0] == 18 && coor[1] == 5) return 63;
else if (coor[0] == 0 && coor[1] == 11) return 45;
else if (coor[0] == 0 && coor[1] == 7) return 15;
else if (coor[0] == 0 && coor[1] == 2) return 2;
else if (coor[0] == 0 && coor[1] == 17) return 58;
else if (coor[0] == 0 && coor[1] == 1) return 0;
else if (coor[0] == 0 && coor[1] == 18) return 62;
else if (coor[0] == 0 && coor[1] == 0) return 0;
else if (coor[0] == 0 && coor[1] == 16) return 47;
else if (coor[0] == 0 && coor[1] == 13) return 30;
else if (coor[0] == 0 && coor[1] == 6) return 10;
else if (coor[0] == 0 && coor[1] == 3) return 4;
else if (coor[0] == 0 && coor[1] == 9) return 28;
else if (coor[0] == 0 && coor[1] == 12) return 55;
else if (coor[0] == 0 && coor[1] == 15) return 61;
else if (coor[0] == 0 && coor[1] == 14) return 47;
else if (coor[0] == 0 && coor[1] == 8) return 21;
else if (coor[0] == 0 && coor[1] == 4) return 6;
else if (coor[0] == 0 && coor[1] == 10) return 36;
else if (coor[0] == 0 && coor[1] == 19) return 32;
else if (coor[0] == 0 && coor[1] == 5) return 8;
else if (coor[0] == 13 && coor[1] == 11) return 15;
else if (coor[0] == 13 && coor[1] == 7) return 36;
else if (coor[0] == 13 && coor[1] == 17) return 27;
else if (coor[0] == 13 && coor[1] == 2) return 56;
else if (coor[0] == 13 && coor[1] == 1) return 38;
else if (coor[0] == 13 && coor[1] == 18) return 43;
else if (coor[0] == 13 && coor[1] == 0) return 31;
else if (coor[0] == 13 && coor[1] == 16) return 39;
else if (coor[0] == 13 && coor[1] == 13) return 22;
else if (coor[0] == 13 && coor[1] == 6) return 45;
else if (coor[0] == 13 && coor[1] == 3) return 35;
else if (coor[0] == 13 && coor[1] == 9) return 20;
else if (coor[0] == 13 && coor[1] == 12) return 54;
else if (coor[0] == 13 && coor[1] == 15) return 51;
else if (coor[0] == 13 && coor[1] == 14) return 17;
else if (coor[0] == 13 && coor[1] == 8) return 63;
else if (coor[0] == 13 && coor[1] == 4) return 42;
else if (coor[0] == 13 && coor[1] == 10) return 54;
else if (coor[0] == 13 && coor[1] == 19) return 23;
else if (coor[0] == 13 && coor[1] == 5) return 29;
else if (coor[0] == 16 && coor[1] == 11) return 8;
else if (coor[0] == 16 && coor[1] == 7) return 8;
else if (coor[0] == 16 && coor[1] == 2) return 57;
else if (coor[0] == 16 && coor[1] == 17) return 4;
else if (coor[0] == 16 && coor[1] == 1) return 48;
else if (coor[0] == 16 && coor[1] == 18) return 29;
else if (coor[0] == 16 && coor[1] == 0) return 48;
else if (coor[0] == 16 && coor[1] == 16) return 49;
else if (coor[0] == 16 && coor[1] == 13) return 8;
else if (coor[0] == 16 && coor[1] == 6) return 53;
else if (coor[0] == 16 && coor[1] == 3) return 38;
else if (coor[0] == 16 && coor[1] == 9) return 14;
else if (coor[0] == 16 && coor[1] == 12) return 58;
else if (coor[0] == 16 && coor[1] == 14) return 53;
else if (coor[0] == 16 && coor[1] == 15) return 25;
else if (coor[0] == 16 && coor[1] == 8) return 18;
else if (coor[0] == 16 && coor[1] == 4) return 40;
else if (coor[0] == 16 && coor[1] == 10) return 18;
else if (coor[0] == 16 && coor[1] == 19) return 43;
else if (coor[0] == 16 && coor[1] == 5) return 44;
else if (coor[0] == 6 && coor[1] == 11) return 63;
else if (coor[0] == 6 && coor[1] == 7) return 49;
else if (coor[0] == 6 && coor[1] == 17) return 48;
else if (coor[0] == 6 && coor[1] == 2) return 24;
else if (coor[0] == 6 && coor[1] == 1) return 16;
else if (coor[0] == 6 && coor[1] == 18) return 6;
else if (coor[0] == 6 && coor[1] == 0) return 11;
else if (coor[0] == 6 && coor[1] == 13) return 63;
else if (coor[0] == 6 && coor[1] == 16) return 52;
else if (coor[0] == 6 && coor[1] == 6) return 32;
else if (coor[0] == 6 && coor[1] == 3) return 33;
else if (coor[0] == 6 && coor[1] == 9) return 50;
else if (coor[0] == 6 && coor[1] == 12) return 51;
else if (coor[0] == 6 && coor[1] == 14) return 12;
else if (coor[0] == 6 && coor[1] == 15) return 21;
else if (coor[0] == 6 && coor[1] == 8) return 63;
else if (coor[0] == 6 && coor[1] == 4) return 43;
else if (coor[0] == 6 && coor[1] == 10) return 61;
else if (coor[0] == 6 && coor[1] == 19) return 16;
else if (coor[0] == 6 && coor[1] == 5) return 54;
else if (coor[0] == 3 && coor[1] == 11) return 53;
else if (coor[0] == 3 && coor[1] == 7) return 40;
else if (coor[0] == 3 && coor[1] == 17) return 23;
else if (coor[0] == 3 && coor[1] == 2) return 9;
else if (coor[0] == 3 && coor[1] == 1) return 5;
else if (coor[0] == 3 && coor[1] == 18) return 60;
else if (coor[0] == 3 && coor[1] == 0) return 5;
else if (coor[0] == 3 && coor[1] == 13) return 34;
else if (coor[0] == 3 && coor[1] == 16) return 37;
else if (coor[0] == 3 && coor[1] == 6) return 32;
else if (coor[0] == 3 && coor[1] == 3) return 14;
else if (coor[0] == 3 && coor[1] == 9) return 59;
else if (coor[0] == 3 && coor[1] == 12) return 59;
else if (coor[0] == 3 && coor[1] == 15) return 38;
else if (coor[0] == 3 && coor[1] == 14) return 49;
else if (coor[0] == 3 && coor[1] == 8) return 49;
else if (coor[0] == 3 && coor[1] == 4) return 19;
else if (coor[0] == 3 && coor[1] == 10) return 39;
else if (coor[0] == 3 && coor[1] == 19) return 40;
else if (coor[0] == 3 && coor[1] == 5) return 25;
else if (coor[0] == 9 && coor[1] == 11) return 27;
else if (coor[0] == 9 && coor[1] == 7) return 51;
else if (coor[0] == 9 && coor[1] == 17) return 56;
else if (coor[0] == 9 && coor[1] == 2) return 48;
else if (coor[0] == 9 && coor[1] == 1) return 37;
else if (coor[0] == 9 && coor[1] == 18) return 37;
else if (coor[0] == 9 && coor[1] == 0) return 29;
else if (coor[0] == 9 && coor[1] == 16) return 11;
else if (coor[0] == 9 && coor[1] == 13) return 19;
else if (coor[0] == 9 && coor[1] == 6) return 51;
else if (coor[0] == 9 && coor[1] == 3) return 60;
else if (coor[0] == 9 && coor[1] == 9) return 60;
else if (coor[0] == 9 && coor[1] == 12) return 41;
else if (coor[0] == 9 && coor[1] == 14) return 23;
else if (coor[0] == 9 && coor[1] == 15) return 51;
else if (coor[0] == 9 && coor[1] == 8) return 37;
else if (coor[0] == 9 && coor[1] == 4) return 42;
else if (coor[0] == 9 && coor[1] == 10) return 63;
else if (coor[0] == 9 && coor[1] == 19) return 31;
else if (coor[0] == 9 && coor[1] == 5) return 60;
else if (coor[0] == 12 && coor[1] == 11) return 4;
else if (coor[0] == 12 && coor[1] == 7) return 57;
else if (coor[0] == 12 && coor[1] == 17) return 47;
else if (coor[0] == 12 && coor[1] == 2) return 52;
else if (coor[0] == 12 && coor[1] == 1) return 34;
else if (coor[0] == 12 && coor[1] == 18) return 8;
else if (coor[0] == 12 && coor[1] == 0) return 56;
else if (coor[0] == 12 && coor[1] == 16) return 18;
else if (coor[0] == 12 && coor[1] == 13) return 51;
else if (coor[0] == 12 && coor[1] == 6) return 54;
else if (coor[0] == 12 && coor[1] == 3) return 60;
else if (coor[0] == 12 && coor[1] == 9) return 42;
else if (coor[0] == 12 && coor[1] == 12) return 17;
else if (coor[0] == 12 && coor[1] == 15) return 5;
else if (coor[0] == 12 && coor[1] == 14) return 39;
else if (coor[0] == 12 && coor[1] == 8) return 53;
else if (coor[0] == 12 && coor[1] == 4) return 42;
else if (coor[0] == 12 && coor[1] == 10) return 7;
else if (coor[0] == 12 && coor[1] == 19) return 60;
else if (coor[0] == 12 && coor[1] == 5) return 61;
else if (coor[0] == 15 && coor[1] == 11) return 51;
else if (coor[0] == 15 && coor[1] == 7) return 45;
else if (coor[0] == 15 && coor[1] == 17) return 62;
else if (coor[0] == 15 && coor[1] == 2) return 56;
else if (coor[0] == 15 && coor[1] == 1) return 56;
else if (coor[0] == 15 && coor[1] == 18) return 18;
else if (coor[0] == 15 && coor[1] == 0) return 62;
else if (coor[0] == 15 && coor[1] == 13) return 54;
else if (coor[0] == 15 && coor[1] == 16) return 24;
else if (coor[0] == 15 && coor[1] == 6) return 35;
else if (coor[0] == 15 && coor[1] == 3) return 50;
else if (coor[0] == 15 && coor[1] == 9) return 54;
else if (coor[0] == 15 && coor[1] == 12) return 7;
else if (coor[0] == 15 && coor[1] == 14) return 60;
else if (coor[0] == 15 && coor[1] == 15) return 0;
else if (coor[0] == 15 && coor[1] == 8) return 20;
else if (coor[0] == 15 && coor[1] == 4) return 28;
else if (coor[0] == 15 && coor[1] == 19) return 36;
else if (coor[0] == 15 && coor[1] == 10) return 13;
else if (coor[0] == 15 && coor[1] == 5) return 55;
else if (coor[0] == 14 && coor[1] == 11) return 17;
else if (coor[0] == 14 && coor[1] == 7) return 58;
else if (coor[0] == 14 && coor[1] == 2) return 31;
else if (coor[0] == 14 && coor[1] == 17) return 23;
else if (coor[0] == 14 && coor[1] == 1) return 51;
else if (coor[0] == 14 && coor[1] == 18) return 6;
else if (coor[0] == 14 && coor[1] == 0) return 48;
else if (coor[0] == 14 && coor[1] == 13) return 18;
else if (coor[0] == 14 && coor[1] == 16) return 52;
else if (coor[0] == 14 && coor[1] == 6) return 13;
else if (coor[0] == 14 && coor[1] == 3) return 55;
else if (coor[0] == 14 && coor[1] == 9) return 58;
else if (coor[0] == 14 && coor[1] == 12) return 8;
else if (coor[0] == 14 && coor[1] == 15) return 59;
else if (coor[0] == 14 && coor[1] == 14) return 39;
else if (coor[0] == 14 && coor[1] == 8) return 46;
else if (coor[0] == 14 && coor[1] == 4) return 21;
else if (coor[0] == 14 && coor[1] == 10) return 8;
else if (coor[0] == 14 && coor[1] == 19) return 62;
else if (coor[0] == 14 && coor[1] == 5) return 48;
else if (coor[0] == 8 && coor[1] == 11) return 34;
else if (coor[0] == 8 && coor[1] == 7) return 29;
else if (coor[0] == 8 && coor[1] == 2) return 39;
else if (coor[0] == 8 && coor[1] == 17) return 33;
else if (coor[0] == 8 && coor[1] == 1) return 29;
else if (coor[0] == 8 && coor[1] == 18) return 14;
else if (coor[0] == 8 && coor[1] == 0) return 22;
else if (coor[0] == 8 && coor[1] == 16) return 17;
else if (coor[0] == 8 && coor[1] == 13) return 48;
else if (coor[0] == 8 && coor[1] == 6) return 61;
else if (coor[0] == 8 && coor[1] == 3) return 50;
else if (coor[0] == 8 && coor[1] == 9) return 61;
else if (coor[0] == 8 && coor[1] == 12) return 52;
else if (coor[0] == 8 && coor[1] == 15) return 19;
else if (coor[0] == 8 && coor[1] == 14) return 22;
else if (coor[0] == 8 && coor[1] == 8) return 63;
else if (coor[0] == 8 && coor[1] == 4) return 62;
else if (coor[0] == 8 && coor[1] == 19) return 5;
else if (coor[0] == 8 && coor[1] == 10) return 44;
else if (coor[0] == 8 && coor[1] == 5) return 44;
else if (coor[0] == 4 && coor[1] == 11) return 62;
else if (coor[0] == 4 && coor[1] == 7) return 51;
else if (coor[0] == 4 && coor[1] == 17) return 43;
else if (coor[0] == 4 && coor[1] == 2) return 13;
else if (coor[0] == 4 && coor[1] == 1) return 7;
else if (coor[0] == 4 && coor[1] == 18) return 57;
else if (coor[0] == 4 && coor[1] == 0) return 7;
else if (coor[0] == 4 && coor[1] == 16) return 27;
else if (coor[0] == 4 && coor[1] == 13) return 41;
else if (coor[0] == 4 && coor[1] == 6) return 42;
else if (coor[0] == 4 && coor[1] == 3) return 20;
else if (coor[0] == 4 && coor[1] == 9) return 41;
else if (coor[0] == 4 && coor[1] == 12) return 41;
else if (coor[0] == 4 && coor[1] == 15) return 23;
else if (coor[0] == 4 && coor[1] == 14) return 20;
else if (coor[0] == 4 && coor[1] == 8) return 61;
else if (coor[0] == 4 && coor[1] == 4) return 27;
else if (coor[0] == 4 && coor[1] == 19) return 21;
else if (coor[0] == 4 && coor[1] == 10) return 57;
else if (coor[0] == 4 && coor[1] == 5) return 34;
else if (coor[0] == 19 && coor[1] == 11) return 56;
else if (coor[0] == 19 && coor[1] == 7) return 0;
else if (coor[0] == 19 && coor[1] == 17) return 20;
else if (coor[0] == 19 && coor[1] == 2) return 39;
else if (coor[0] == 19 && coor[1] == 1) return 62;
else if (coor[0] == 19 && coor[1] == 18) return 5;
else if (coor[0] == 19 && coor[1] == 0) return 33;
else if (coor[0] == 19 && coor[1] == 16) return 15;
else if (coor[0] == 19 && coor[1] == 13) return 0;
else if (coor[0] == 19 && coor[1] == 6) return 19;
else if (coor[0] == 19 && coor[1] == 3) return 43;
else if (coor[0] == 19 && coor[1] == 9) return 32;
else if (coor[0] == 19 && coor[1] == 12) return 61;
else if (coor[0] == 19 && coor[1] == 15) return 39;
else if (coor[0] == 19 && coor[1] == 14) return 3;
else if (coor[0] == 19 && coor[1] == 8) return 7;
else if (coor[0] == 19 && coor[1] == 4) return 35;
else if (coor[0] == 19 && coor[1] == 10) return 12;
else if (coor[0] == 19 && coor[1] == 19) return 35;
else if (coor[0] == 19 && coor[1] == 5) return 48;
else if (coor[0] == 10 && coor[1] == 11) return 36;
else if (coor[0] == 10 && coor[1] == 7) return 41;
else if (coor[0] == 10 && coor[1] == 17) return 51;
else if (coor[0] == 10 && coor[1] == 2) return 58;
else if (coor[0] == 10 && coor[1] == 1) return 46;
else if (coor[0] == 10 && coor[1] == 18) return 39;
else if (coor[0] == 10 && coor[1] == 0) return 37;
else if (coor[0] == 10 && coor[1] == 13) return 51;
else if (coor[0] == 10 && coor[1] == 16) return 17;
else if (coor[0] == 10 && coor[1] == 6) return 62;
else if (coor[0] == 10 && coor[1] == 3) return 40;
else if (coor[0] == 10 && coor[1] == 9) return 26;
else if (coor[0] == 10 && coor[1] == 12) return 23;
else if (coor[0] == 10 && coor[1] == 15) return 12;
else if (coor[0] == 10 && coor[1] == 14) return 14;
else if (coor[0] == 10 && coor[1] == 8) return 61;
else if (coor[0] == 10 && coor[1] == 4) return 58;
else if (coor[0] == 10 && coor[1] == 19) return 27;
else if (coor[0] == 10 && coor[1] == 10) return 42;
else if (coor[0] == 10 && coor[1] == 5) return 45;
else if (coor[0] == 5 && coor[1] == 11) return 52;
else if (coor[0] == 5 && coor[1] == 7) return 63;
else if (coor[0] == 5 && coor[1] == 2) return 18;
else if (coor[0] == 5 && coor[1] == 17) return 38;
else if (coor[0] == 5 && coor[1] == 1) return 11;
else if (coor[0] == 5 && coor[1] == 18) return 58;
else if (coor[0] == 5 && coor[1] == 0) return 9;
else if (coor[0] == 5 && coor[1] == 13) return 28;
else if (coor[0] == 5 && coor[1] == 16) return 57;
else if (coor[0] == 5 && coor[1] == 6) return 53;
else if (coor[0] == 5 && coor[1] == 3) return 26;
else if (coor[0] == 5 && coor[1] == 9) return 59;
else if (coor[0] == 5 && coor[1] == 12) return 57;
else if (coor[0] == 5 && coor[1] == 15) return 46;
else if (coor[0] == 5 && coor[1] == 14) return 47;
else if (coor[0] == 5 && coor[1] == 8) return 43;
else if (coor[0] == 5 && coor[1] == 4) return 35;
else if (coor[0] == 5 && coor[1] == 19) return 29;
else if (coor[0] == 5 && coor[1] == 10) return 36;
else if (coor[0] == 5 && coor[1] == 5) return 44;


    // int tileYIndex = coor[1]  / peCols;
//     int XwithinPEtile = (coor[0] + tileYIndex * peRotate) % peRows;
//     int YwithinPEtile = coor[1] % (peCols / peStride);
//     int subtileY = (coor[1] % peCols) / (peCols / peStride);
//     int peNum = XwithinPEtile * peStride + YwithinPEtile * peStride * peRows + subtileY;
//     CkAssert(peNum < CkNumPes());
//     return peNum;
  }

  int procNum(int arrayHdl, const CkArrayIndex &idx) {
    const int *coor = idx.data();
    return map(coor);
  }

  int pesInPanel(CkIndex2D index) {
    if (index.x < index.y)
      return peRows;
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
      int homePE = map(coor);
      pesInRow.erase(homePE);

      return pesInRow.size();
    }
  }

private:
  const int peRows, peCols, peRotate, peStride, numBlks;
};

/// Mapping for BlockScheduler - created as an array so group section multicast
/// can be used
struct OnePerPE : public CkArrayMap {
  OnePerPE() { }
  //int registerArray(CkArrayIndexMax& numElements,CkArrayID aid) {
  //CkAssert(CkNumPes() == numElements.index[0]);
  //return 0;
  //}
  int procNum(int arrayHdl, const CkArrayIndex &elt) {
    assert(*(int*)elt.data() < CkNumPes() &&
	   *(int*)elt.data() >= 0);
    //CkPrintf("index = %d\n", *(int*)elt.data());
    //fflush(stdout);
    return *(int*)elt.data();
  }
};
#endif
