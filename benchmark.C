#include "luConfig.h"
#include "benchmark.decl.h"
#include <algorithm>

class Benchmark : public CBase_Benchmark {
  LUConfig luCfg;
  int runs;

public:
  Benchmark(CkArgMsg* m) {
    if (m->argc<4) {
      CkPrintf("Usage: %s <matrix size> <block size> <mem threshold> [<pivot batch size> <mapping scheme> [<peTileRows> <peTileCols>] ]\n", m->argv[0]);
      CkExit();
    }

    /// Parse the command line and accept user input
    luCfg.matrixSize     = atoi(m->argv[1]);
    luCfg.blockSize      = atoi(m->argv[2]);
    luCfg.memThreshold   = atoi(m->argv[3]);

    if (m->argc >= 5)
      luCfg.pivotBatchSize = atoi( m->argv[4] );
    else
      luCfg.pivotBatchSize = luCfg.blockSize / 4;

    if (m->argc >= 6) {
      luCfg.mappingScheme = atoi(m->argv[5]);
      if (luCfg.mappingScheme == 3) {
        if (m->argc >= 8) {
          luCfg.peTileRows = atoi( m->argv[6] );
          luCfg.peTileCols = atoi( m->argv[7] );
        }
        if (m->argc < 9)
          luCfg.peTileRotate = 0;
        else
          luCfg.peTileRotate = atoi(m->argv[8]);
        int peTileSize = luCfg.peTileRows * luCfg.peTileCols;
        if (peTileSize > CkNumPes())
          CkAbort("The PE tile dimensions are too big for the num of PEs available!");
        if (peTileSize < CkNumPes())
          CkPrintf("WARNING: Configured to use a PE tile size(%dx%d)"
                   "(for 2D tile mapping) that does not use all the PEs(%d)\n",
                   luCfg.peTileRows, luCfg.peTileCols, CkNumPes() );
      }
    }
    else
      luCfg.mappingScheme = 1;

    if (luCfg.matrixSize % luCfg.blockSize!=0) {
      CkPrintf("The matrix size %d should be a multiple of block size %d!\n",
               luCfg.matrixSize, luCfg.blockSize);
      CkExit();
    }

    luCfg.numBlocks = luCfg.matrixSize / luCfg.blockSize;

    runs = 1;
    for (int i = 0; i < runs; i++) {
      ckout << "Test harness instigating solver" << endl;
      CkCallback cb(CkIndex_Benchmark::finished(), thisProxy);
      CProxy_LUSolver solver = CProxy_LUSolver::ckNew(luCfg, cb);
    }
  }

  void finished() {
    ckout << "Test harness received callback" << endl;
    if (--runs == 0)
      CkExit();
  }
};

#include "benchmark.def.h"
