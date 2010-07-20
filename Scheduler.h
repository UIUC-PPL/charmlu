#ifndef SCHEDULER_H
#define SCHEDULER_H

#include <vector>
#include <list>
#include <vector>
#include <climits>
#include "pup_stl.h"
#include "c_common.h"
#include "scheduler.decl.h"

using namespace std;

struct JMessage {
  CProxyElement_LUBlk block;
  int x, y;
  int cx, cy;
  double *first;
  double *second;
  double *LU;

  int fsize, ssize;
  int step;

  int getSize();
};

class Scheduler : public CBase_Scheduler {
public:
  Scheduler() : count(0), total_size(0), sch_count(0),
		msg_count(0), gpu_count(0), cpu_count(0),
		GPUworking(false), cpuStatus(CkNumPes()), gpu_state(0),
                maxRow(0), maxCol(0) {}
  Scheduler(CkMigrateMessage* msg) {}

  list<JMessage*> mapMsg;
  vector<bool> cpuStatus;
  map<int, list<JMessage*> > rowMap, colMap;
  int maxRow, maxCol;
  int maxRowSize, maxColSize;

  int count, total_size, sch_count, msg_count, gpu_count, cpu_count;
  int gpu_state;
  bool GPUworking;

  void tryAgain(int a);
  void cpuFree(int cpu);
  
  void finishedGPU(list<JMessage> msgs);

  void haveWork(CProxyElement_LUBlk block, double** first, 
		double** second, double** LU, int x, int y, int cx,
		int cy, int fsize, int ssize, int step);
  
  void checkIn();

  void findNewMax();
  list<JMessage*>* findLargestAgglom();
  void removeItemXY(list<JMessage*>& list1, JMessage* &msg);

  // Just for testing
  void square_dgemm(int N, double *A, double *B, double *C);

private:
  int findLargestInMap(map<int, list<JMessage*> >);
};

#endif
