#include "hello.decl.h"
#include <stdio.h>
#include <vector>
#include <algorithm>

/*readonly*/ CProxy_Main mainProxy;

/*mainchare*/
class Main : public CBase_Main
{
public:
  Main(CkArgMsg* m)
  {
    //Process command-line arguments
    if (m->argc < 2)
        CkAbort("Usage: <exe> blockSize\nwhere, blocksize is the size of 2D array blocks that should be allocated\n");
    blocksize = atoi(m->argv[1]);
    delete m;

    // Other initialization
    numFinished = 0;
    numblocks.resize(CkNumPes(), 0);

    //Start the computation
    CkPrintf("Running memory fragmentation test on %d processors\n", CkNumPes());
    mainProxy = thisProxy;

    CProxy_Hello grp = CProxy_Hello::ckNew();
    grp.fillAndTestMem(blocksize);
  }

  void done(int pe, int numAllocated)
  {
    numblocks[pe] = numAllocated;

    if (++numFinished == CkNumPes()) {

      CkPrintf("\nblock size: %d x %d doubles", blocksize, blocksize);
      CkPrintf("\nmem size  : %d KB", 8.0 * blocksize*blocksize / 1024);

      CkPrintf("\nmin number of blocks that fit on any given PE: %d",
                    * std::min_element(numblocks.begin(), numblocks.end()) );
      CkPrintf("\nmax number of blocks that fit on any given PE: %d",
                    * std::max_element(numblocks.begin(), numblocks.end()) );

      for (int i=0; i<numblocks.size(); i++)
        CkPrintf("%d: %d\n",i, numblocks[i]);

      CkExit();
    }
  }

private:
  int blocksize;
  std::vector<int> numblocks;
  int numFinished;
};


class Hello : public CBase_Hello
{
public:
  void fillAndTestMem(int blocksize)
  {
    int arrsize = blocksize * blocksize;
    std::vector<double*> blocks;
    int numAllocated;
    try {
      while(1) { blocks.push_back(new double[arrsize]); }
    }
    catch (...) {
      numAllocated = blocks.size();
      for (int i=0; i< blocks.size(); i++)
        delete [] blocks[i];
      blocks.clear();
    }

    mainProxy.done(CkMyPe(), numAllocated);
  }
};

#include "hello.def.h"
