#include "benchmark.decl.h"
#include <vector>
#include <algorithm>
#include "acml.h"

#define BLAS_NOTRANSPOSE 'N'

CProxy_Main mainProxy;
int blockSize;
int numIter;

struct Main : public CBase_Main {
  int count;
  CProxy_dgemmTest dgemms;
  std::vector<double> peTimes;

  Main(CkArgMsg *m) : count(0), peTimes(CkNumPes()) {
    if (m->argc < 3) {
      CkPrintf("usage: blockSize numIter\n");
      CkExit();
    }

    blockSize = atoi(m->argv[1]);
    numIter = atoi(m->argv[2]);

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

  double *m1, *m2, *m3;

  dgemmTest() {
    m1 = new double[blockSize*blockSize];
    m2 = new double[blockSize*blockSize];
    m3 = new double[blockSize*blockSize];
  }

  void testAndSync() {

    double startTestTime = CmiWallTimer();

    for(int i=0; i < numIter; i++) {
        dgemm(BLAS_NOTRANSPOSE, BLAS_NOTRANSPOSE,
            blockSize, blockSize, blockSize,
            -1.0, m1,
            blockSize, m2, blockSize,
            1.0, m3, blockSize);
    }

    double totalTestTime = CmiWallTimer() - startTestTime;

    delete[] m1;
    delete[] m2;
    delete[] m3;

    mainProxy.finishedTest(CkMyPe(), totalTestTime);
  }
};

#include "benchmark.def.h"
