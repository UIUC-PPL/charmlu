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


/** do a space filling curve style allocation from the bottom right to the upper left. */
class LUSnakeMap: public LUMap {
public:
  int numBlks, BLKSIZE;
  LUSnakeMap(int _numBlks, int _BLKSIZE) : numBlks(_numBlks), BLKSIZE(_BLKSIZE) {}
  int procNum(int arrayHdl, const CkArrayIndex &idx) {
    int *coor = (int *)idx.data();
    int x = coor[0];
    int y = coor[1];

    int numProcs = CkNumPes();

    const int numsteps=numBlks;
   
    int p=0;
    for(int step=numsteps-1;step>=0;step--){
      
      // go along row
      for(int i=1;i<numsteps-step;i++){
	int curr_x = i + step;
	int curr_y = step;
	if(x==curr_x && y==curr_y)
	  return p % numProcs;
	p++;
      }
      
      // visit corner
      if(x==step && y==step)
	return p % numProcs;
      p++;
      
      // go along column
      for(int i=1;i<numsteps-step;i++){
	int curr_y = i + step;
	int curr_x = step;
	if(x==curr_x && y==curr_y)
	  return p % numProcs;
	p++;
      }
    }
    
    CkAbort("Mapping code has a bug in it\n");
    return -1;
  }
};

/** do an allocation that results in almost identical numbers of trailing updates per PE */
class LUBalancedSnakeMap: public LUMap {
public:
  
  int mappingSize;
  int *mapping;
  int *peLoads;
  int numBlks, BLKSIZE;
  

  void setMapping(int x, int y, int pe){
    CkAssert(x*numBlks+y < mappingSize);
    mapping[x*numBlks+y] = pe;
    peLoads[pe] += workLoad(x, y);
  }

  int getMapping(int x, int y){
    CkAssert(x*numBlks+y < mappingSize);
    return mapping[x*numBlks+y];
  }

  /** build and store the mapping once */
  LUBalancedSnakeMap(int _numBlks, int _BLKSIZE) : numBlks(_numBlks), BLKSIZE(_BLKSIZE) {
    int numProcs = CkNumPes();

    mappingSize = numBlks*numBlks;
    mapping = new int[mappingSize];    

    peLoads = new int[numProcs];
    for(int i=0;i<numProcs;i++)
      peLoads[i]=0;
    
    for(int step = numBlks-1; step >= 0; step--) {
        // visit corner & column
        for(int x = step; x < numBlks; x++)
            setMapping(x, step, minLoadedPE() );
        // go along row
        for(int y = step+1; y < numBlks; y++)
            setMapping(step, y, minLoadedPE() );
    }
  }
  
  int minLoadedPE(){
    int minLoadFound = 1000000000;
    int minPEFound = CkNumPes()-1;
    for(int p=CkNumPes()-1; p>=0; p--){
      if(peLoads[p] < minLoadFound) {
	minPEFound = p;
	minLoadFound = peLoads[p];
      }
    }
    return minPEFound;	  
  }
  
  int workLoad(int x, int y){
    if(x<y)
      return x+1;
    else 
      return y+1;
  }
  
  int procNum(int arrayHdl, const CkArrayIndex &idx) {
    int *coor = (int *)idx.data();
    int x = coor[0];
    int y = coor[1];
    return getMapping(x,y);
  }

};

/** do an allocation that results in almost identical numbers of trailing updates per PE */
class LUBalancedSnakeMap2: public LUMap {
public:
  int mappingSize;
  int *mapping;
  int *peLoads;
  int stateN;
  int numBlks, BLKSIZE;

  void setMapping(int x, int y, int pe){
    CkAssert(y*numBlks+x < mappingSize);
    mapping[y*numBlks+x] = pe;
    peLoads[pe] += workLoad(x, y);
  }

  int getMapping(int x, int y){
    CkAssert(y*numBlks+x < mappingSize);
    return mapping[y*numBlks+x];
  }

  /** build and store the mapping once */
  LUBalancedSnakeMap2(int _numBlks, int _BLKSIZE) : numBlks(_numBlks), BLKSIZE(_BLKSIZE) {
    stateN = 0;
    int numProcs = CkNumPes();

    mappingSize = numBlks*numBlks;
    mapping = new int[mappingSize];  

    peLoads = new int[numProcs];
    for (int i = 0; i < numProcs; i++)
      peLoads[i]=0;
    
    const int numsteps = numBlks;

    for(int step = numsteps - 1; step >= 2; step--){
      // go along row
      for(int i = 1; i < numsteps - step; i++){
	int x = i + step;
	int y = step;
	int minLoaded = minLoadedPE();
	setMapping(x, y, minLoaded);
      }
      
      // visit corner & column
      for(int i = 0; i < numsteps - step; i++){
	int y = i + step;
	int x = step;
	int minLoaded = minLoadedPE();
	setMapping(x, y, minLoaded);
      }
    }

    // go along first two rows
    for (int x = 1; x < numsteps; x++){
      int minLoaded = minLoadedPE();
      setMapping(x, 0, minLoaded);
      setMapping(x, 1, minLoaded);
    }
    
    // visit first corner & first two columns
    for (int y = 0; y < numsteps; y++) {
      int minLoaded = minLoadedPE();
      setMapping(0, y, minLoaded);
      if (y != 0)
	setMapping(1, y, minLoaded);
    }

  }
  
  int minLoadedPE() {
    int minLoadFound = 1000000000;
    int minPEFound = CkNumPes()-1;
    for(int p = CkNumPes() - 1; p >= 0; p--){
      if(peLoads[p] < minLoadFound) {
	minPEFound = p;
	minLoadFound = peLoads[p];
      }
    }

    if (stateN == minPEFound) {
      int proc = stateN++ % CkNumPes();
      stateN = proc;
      return proc;
    }

    stateN = minPEFound;
    return minPEFound;	  
  }
  
  int workLoad(int x, int y){
    if (x < y)
      return x+1;
    else 
      return y+1;
  }
  
  int procNum(int arrayHdl, const CkArrayIndex &idx) {
    int *coor = (int *)idx.data();
    int x = coor[0];
    int y = coor[1];
    return getMapping(x,y);
  }

};

class BlockCyclicMap: public LUMap {
public:
  BlockCyclicMap() {}
  int procNum(int arrayHdl, const CkArrayIndex &idx) {
    int *coor = (int *)idx.data();
    int x = coor[1];
    int y = coor[0];

    int numProcs = CkNumPes();
    int numNodes = CkNumNodes();
    //int numNodes = numProcs/4;

    int procPerNode = numProcs/numNodes;
    //int procPerNode = 4;
    CkAssert(numProcs%numNodes==0);

    //assume a procPerNode X numNodes grid of processors

    int penum = (x%procPerNode)*numNodes+(y%numNodes);
    return penum;
  }
};

class RealBlockCyclicMap : public LUMap {
    int r, num_blocks;
public:
    RealBlockCyclicMap(int r_, int num_blocks_) : r(r_), num_blocks(num_blocks_) {}
    int procNum(int arrayHdl, const CkArrayIndex &idx) {
	int *coor = (int *)idx.data();
	int m = coor[1]*num_blocks + coor[0];
	return (m % (r*CkNumPes()))/r;
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
    int tileYIndex = coor[1]  / peCols;
    int XwithinPEtile = (coor[0] + tileYIndex * peRotate) % peRows;
    int YwithinPEtile = coor[1] % peStride;
    int subtileY = (coor[1] % peCols) / peStride;
    int peNum = XwithinPEtile * peStride + YwithinPEtile * peRows + subtileY;
    CkAssert(peNum < CkNumPes());
    return peNum;
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

#include <utility>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

bool mapValueSort (const std::pair<int, int> &p1, const std::pair<int, int> &p2) {
  return p1.second < p2.second;
}

class StrongScaling1 : public LUMap {
private:
  int numBlks;
  std::map<std::pair<int, int>, int> peBlock;

  int nextPE(int pe) {
    return (pe + 1) % CkNumPes();
  }

public:
    StrongScaling1(int numBlks_) : numBlks(numBlks_) {
      int numPEs = CkNumPes();
      std::vector<int> peWork(numPEs);

      int currentPE = 0;

      for (int i = 0; i < numPEs; i++) {
        peWork[i] = 0;
      }

      for (int x = 0; x < numBlks; x++) {
        for (int y = 0; y < numBlks; y++) {
          peBlock[std::make_pair(x, y)] = -1;
        }
      }

      for (int y = 0; y < numBlks; y++) {
        for (int x = y; x < numBlks; x++) {
          int pe = nextPE(currentPE);
          if (x == y) {
            // Diagonal "A" block
            peWork[pe] += 15 * y;
          } else if (x > y) {
            // Below diagonal "B" block
            peWork[pe] += 10 * y;
          }
          peBlock[std::make_pair(x, y)] = pe;
          currentPE = pe;
        }
        //startPE = (numBlks-y) % numPEs;
      }

      for (int y = 0; y < numBlks; y++) {
        sort(peWork.begin(), peWork.end());
        for (int x = 0; x < y; x++) {
          int pe = nextPE(currentPE);
          peWork[pe] += 5 * x;
          peBlock[std::make_pair(x, y)] = pe;
          currentPE = pe;
        }
      }

      /* ckout << "--" << endl; */
      /* for (int x = 0; x < numBlks; x++) { */
      /*   for (int y = 0; y < numBlks; y++) { */
      /*     ckout << peBlock[std::make_pair(x, y)] << " "; */
      /*   } */
      /*   ckout << endl; */
      /* } */
      /* ckout << "--" << endl; */
    }

    int procNum(int arrayHdl, const CkArrayIndex &idx) {
	int *coor = (int *)idx.data();
        int x = coor[0], y = coor[1];
        int pe = peBlock[std::make_pair(x, y)];
        CkAssert(pe >= 0);
        return pe;
    }

    std::string desc() { return "strong scaling"; }
};

/// Mapping for BlockScheduler - created as an array to hack around
/// inability to do group section multicasts
struct OnePerPE : public CBase_OnePerPE {
  OnePerPE() { }
  int registerArray(CkArrayIndexMax& numElements,CkArrayID aid) {
    CkAssert(CkNumPes() == numElements.index[0]);
    return 0;
  }
  int procNum(int arrayHdl, const CkArrayIndex &elt) {
    return *(int*)elt.data();
  }
};

#endif
