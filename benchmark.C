#include "luConfig.h"
#include "driver.decl.h"
#include <algorithm>
#include <utility>
#include <cmath>

#define _QUOTEIT(x) #x
#define INQUOTES(x) _QUOTEIT(x)

struct Benchmark : public CBase_Benchmark {
  int numIterations;

  Benchmark(CkArgMsg* m)
    : numIterations(1) {
    if (m->argc<4) {
      CkPrintf("Usage: %s <matrix size> <block size> <mem threshold>"
               " [<pivot batch size> <mapping scheme>"
               " [<peTileRows> <peTileCols> <peRotate> <peStride>] ]\n",
               m->argv[0]);
      CkExit();
    }

    LUConfig luCfg;
    /// Parse the command line and accept user input
    luCfg.matrixSize = atoi(m->argv[1]);
    luCfg.blockSize = atoi(m->argv[2]);
    luCfg.memThreshold = atoi(m->argv[3]);

    if (m->argc >= 5)
      luCfg.pivotBatchSize = atoi(m->argv[4]);
    else
      luCfg.pivotBatchSize = luCfg.blockSize / 4;

    if (m->argc >= 6) {
      luCfg.mappingScheme = atoi(m->argv[5]);
      if (luCfg.mappingScheme == 3) {
        if (m->argc < 8) {
          std::pair<int,int> tileDims = computePETileDimensions();
          luCfg.peTileRows = tileDims.first;
          luCfg.peTileCols = tileDims.second;
        }
        else {
          luCfg.peTileRows = atoi(m->argv[6]);
          luCfg.peTileCols = atoi(m->argv[7]);
        }
        if (m->argc < 9)
          luCfg.peTileRotate = 0;
        else
          luCfg.peTileRotate = atoi(m->argv[8]);

        if (m->argc < 10)
          luCfg.peTileStride = luCfg.peTileCols;
        else
          luCfg.peTileStride = atoi(m->argv[9]);

        int peTileSize = luCfg.peTileRows * luCfg.peTileCols;

        if (peTileSize > CkNumPes())
          CkAbort("The PE tile dimensions are too big for the num of PEs available!");

        if (peTileSize < CkNumPes())
          CkPrintf("WARNING: Configured to use a PE tile size (%d x %d)"
                   " (for 2D tile mapping) that does not use all the PEs(%d)\n",
                   luCfg.peTileRows, luCfg.peTileCols, CkNumPes());
      }
    }
    else
      luCfg.mappingScheme = 1;

    if (luCfg.matrixSize >= std::numeric_limits<short>::max() &&
        sizeof(CMK_REFNUM_TYPE) != sizeof(int)) {
      CkPrintf("Refnum size too small for large matrices."
               " Compile charm with build option --with-refnum-type=int\n");
      CkExit();
    }

    if (luCfg.matrixSize % luCfg.blockSize != 0) {
      CkPrintf("The matrix size %d should be a multiple of block size %d!\n",
               luCfg.matrixSize, luCfg.blockSize);
      CkExit();
    }

    luCfg.numBlocks = luCfg.matrixSize / luCfg.blockSize;

    CkPrintf("Running LU compiled from revision: " INQUOTES(LU_REVISION) "\n");
    CkPrintf("Running LU on %d processors (%d nodes): "
             "\n\tMatrix size: %d X %d "
             "\n\tBlock size: %d X %d "
             "\n\tChare Array size: %d X %d"
             "\n\tPivot batch size: %d"
             "\n\tMem Threshold (MB): %d"
             "\n\tSend Limit: %d"
             "\n\tMapping Scheme: %d (%s)\n",
             CkNumPes(), CmiNumNodes(),
             luCfg.matrixSize, luCfg.matrixSize,
             luCfg.blockSize, luCfg.blockSize,
             luCfg.numBlocks, luCfg.numBlocks,
             luCfg.pivotBatchSize,
             luCfg.memThreshold,
             SEND_LIM,
             luCfg.mappingScheme,
             luCfg.mappingScheme == 1 ? "Balanced Snake" :
             (luCfg.mappingScheme==2 ? "Block Cyclic" :
              (luCfg.mappingScheme == 3 ? "2D Tiling" : "Strong Scaling"))
             );
    if (luCfg.mappingScheme == 3)
      CkPrintf("\tMapping PE tile size: %d x %d"
               " rotate %d stride %d \n", luCfg.peTileRows, luCfg.peTileCols,
               luCfg.peTileRotate, luCfg.peTileStride);

#ifdef SCHED_PIVOT_REDN
    CkPrintf("\tPivot Redn Scheduling: On\n");
#else
    CkPrintf("\tPivot Redn Scheduling: Off\n");
#endif
#if defined(CHARMLU_USEG_FROM_BELOW)
    CkPrintf("\tUseg Multicast from: Below\n");
#else
    CkPrintf("\tUseg Multicast from: Diagonal\n");
#endif

    for (int i = 0; i < numIterations; i++) {
      ckout << "Starting solve" << endl;
      CkCallback cb(CkIndex_Benchmark::finished(), thisProxy);
      CProxy_LUSolver solver = CProxy_LUSolver::ckNew(luCfg, cb);
    }
  }

  std::pair<int,int> computePETileDimensions() {
    // Identify two factors that can be used as the tile dimensions for the PE
    // tile
    int factor1 = std::sqrt((double)CkNumPes());
    while (CkNumPes() % factor1 != 0 && factor1 > 0)
      factor1--;
    if (factor1 == 0)
      CkAbort("Couldn't identify a factor of numPEs to"
              " represent the PEs as a 2D tile");

    int factor2 = CkNumPes() / factor1;

    // Set the tile dimensions
    int numPErows = (factor1 >= factor2) ? factor1 : factor2;
    int numPEcols = CkNumPes() / numPErows;

    if (numPErows * numPEcols != CkNumPes())
      CkAbort("The identified tile dimensions dont match the number of PEs!!");

    return std::make_pair(numPErows, numPEcols);
  }

  void finished() {
    ckout << "Finished solve" << endl;
    if (--numIterations == 0)
      CkExit();
  }
};

#include "driver.def.h"
