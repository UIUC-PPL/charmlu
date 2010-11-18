#include <string>

class LUMap : public CkArrayMap {
    virtual std::string desc() { return ""; }
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
    CkAssert(y*numBlks+x < mappingSize);
    mapping[y*numBlks+x] = pe;
    peLoads[pe] += workLoad(x, y);
  }

  int getMapping(int x, int y){
    CkAssert(y*numBlks+x < mappingSize);
    return mapping[y*numBlks+x];
  }

  /** build and store the mapping once */
  LUBalancedSnakeMap(int _numBlks, int _BLKSIZE) : numBlks(_numBlks), BLKSIZE(_BLKSIZE) {
    int numProcs = CkNumPes();

    mappingSize = numBlks*numBlks;
    mapping = new int[mappingSize];    

    peLoads = new int[numProcs];
    for(int i=0;i<numProcs;i++)
      peLoads[i]=0;
    
    const int numsteps=numBlks;
    

    // go along first row
    int step = 0;
    for(int i=1;i<numsteps-step;i++){
      int x = i + step;
      int y = step;
      setMapping(x, y, minLoadedPE() );
    }
    
    // visit first corner & column
    for(int i=0;i<numsteps-step;i++){
      int y = i + step;
      int x = step;
      setMapping(x, y, minLoadedPE() );
    }
    
    for(int step=numsteps-1;step>=1;step--){
      
      // go along row
      for(int i=1;i<numsteps-step;i++){
	int x = i + step;
	int y = step;
	setMapping(x, y, minLoadedPE() );
      }
      
      // visit corner & column
      for(int i=0;i<numsteps-step;i++){
	int y = i + step;
	int x = step;
	setMapping(x, y, minLoadedPE() );
      }
      
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
