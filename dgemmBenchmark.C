#include "benchmark.decl.h"
#include <vector>
#include <algorithm>

#if USE_ACML_H
#include "acml.h"
#define BLAS_NOTRANSPOSE 'N'

#elif USE_ESSL
#define _ESVCPTR
#include <complex>
#include <essl.h>
#define BLAS_NOTRANSPOSE "N"

#endif

CProxy_Main mainProxy;
int blockSize;
int numIter;
int blockNum;

struct Main : public CBase_Main {
  int count;
  CProxy_dgemmTest dgemms;
  std::vector<double> peTimes;

  Main(CkArgMsg *m) : count(0), peTimes(CkNumPes()) {
    if (m->argc < 4) {
      CkPrintf("usage: blockSize numIter blockNum\n");
      CkExit();
    }

    blockSize = atoi(m->argv[1]);
    numIter = atoi(m->argv[2]);
    blockNum = atoi(m->argv[3]);

    CkPrintf("blockSize = %d, numIter = %d\n",blockSize, numIter);

    mainProxy = thisProxy;
    dgemms = CProxy_dgemmTest::ckNew();
    dgemms.testAndSync();
  }

  void finishedTest(int pe, double time) {
    peTimes[pe] = time;
    if (++count == CkNumPes()) {

        for (int i = 0; i < CkNumPes(); ++i) {
            double dgemmDuration    = peTimes[i]/numIter;
            double dgemmFlopCount   = (double)blockSize * (double)blockSize * (double)blockSize * 2.0;
            double dgemmGFlopCount  = dgemmFlopCount / 1000000000.0;
            double dgemmGFlopPerSec = dgemmGFlopCount / dgemmDuration;

            CkPrintf("The dgemm %d x %d takes %g ms and achieves %g GFlop/sec, %g percent of peak on Jaguar\n", blockSize,
            blockSize, dgemmDuration*1000, dgemmGFlopPerSec, dgemmGFlopPerSec/10.3987);
        }

        CkExit();
    }
  }
};

struct dgemmTest : public CBase_dgemmTest {

  double **m1, **m2, **m3;

  dgemmTest() {

    CkPrintf("PE %d using blocknum = %d\n",CkMyPe(),blockNum);
    m1 = new double*[blockNum];
    m2 = new double*[blockNum];
    m3 = new double*[blockNum];    

    for(int i = 0; i < blockNum;i++) {
      m1[i] = new double[blockSize*blockSize];
      m2[i] = new double[blockSize*blockSize];
      m3[i] = new double[blockSize*blockSize];
    
      for(int j = 0; j < blockSize*blockSize; j++) {
       m1[i][j] = 0;
       m2[i][j] = 0;
       m3[i][j] = 0;
       }
    }
  }

  void testAndSync() {

    double startTestTime = CmiWallTimer();

    for(int i=0; i < numIter; i++) {
        int block = i % blockNum;
        dgemm(BLAS_NOTRANSPOSE, BLAS_NOTRANSPOSE,
            blockSize, blockSize, blockSize,
            -1.0, m1[block],
            blockSize, m2[block], blockSize,
            1.0, m3[block], blockSize);
    }

    double totalTestTime = CmiWallTimer() - startTestTime;

    //Don't care about freeing memory
    //delete[] m1;
    //delete[] m2;
    //delete[] m3;

    mainProxy.finishedTest(CkMyPe(), totalTestTime);
  }
};

#include "benchmark.def.h"
