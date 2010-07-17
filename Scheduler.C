/*
 * Scheduler.C
 *
 *  Created on: July 14, 2010
 *      Author: Jonathan
 */

#if USE_CBLAS_H
extern "C" {
#include <cblas.h>
#include <clapack.h>
}

#elif USE_MKL_CBLAS_H
#include "mkl_cblas.h"
#include "mkl_lapack.h"

#elif USE_ACML_H
#include "acml.h"

#elif USE_ACCELERATE_BLAS
#include <Accelerate/Accelerate.h>

#elif USE_ESSL
#define _ESVCPTR
#include <complex>
#include <essl.h>


#else
#error "No BLAS Header files included!"
#endif


#include <list>
#include <vector>
#include <climits>
#include "pup_stl.h"
#include "scheduler.decl.h"
#include "Scheduler.h"


int JMessage::getSize() {
  return fsize + ssize;
}

class CProxy_Cell;
extern CProxy_Main mainProxy;
extern CProxy_GPUWork gpu;
extern CProxy_LUBlk luArrProxy;

using namespace std;

void Scheduler::tryAgain(int a) {
  sch_count++;

  ckout << "calling tryAgain" << endl;

  // GPU offloading
  if (1) {
    int totalSize = 0;
    int numberAgglom = 0;
    list<JMessage*> toOffload;

    if (mapMsg.size() >= MIN_SIZE) {
      while (numberAgglom < MAX_SIZE) {
        list<JMessage*>* msgs = findLargestAgglom();

        if (msgs != NULL) {
          for (list<JMessage*>::iterator iter = msgs->begin();
               iter != msgs->end(); ++iter) {
            int tsize = (*iter)->fsize + (*iter)->ssize;

            if (tsize + totalSize < MAX_SIZE) {
              toOffload.push_back(*iter);

              mapMsg.remove(*iter);
              rowMap[(*iter)->x].remove(*iter);
              colMap[(*iter)->y].remove(*iter);

              numberAgglom++;
            } else {
              break;
            }
          }
        } else {
          break;
        }
      }

      ckout << "numberAgglom = " << numberAgglom << ", queue size = " << mapMsg.size() << endl; 

      if (numberAgglom > 0) {
        // run on GPU
        gpu_count++;
        total_size += toOffload.size();

        int partsize = 0, k = 0;

        list<JMessage> toSend;

        for (list<JMessage*>::iterator iter2 = toOffload.begin(); 
             iter2 != toOffload.end(); ++iter2) {
          toSend.push_back(**iter2);
        }

        GPUworking = true;

        CkEntryOptions opts;
        opts.setPriority(-1000);

        //ckout << "running " << numberAgglom  << " msg on GPU, " << "queue size is " << mapMsg.size() << endl;

        gpu[CkMyPe()].gpu_offload(toSend, &opts);
      }
    }
  }

  // If there is work to do, run one message on CPU
  // for (int i = 0; i < MAX_CPU_SIZE && mapMsg.size() > 0; i++)
  if (mapMsg.size() > 0) {
    cpu_count++;

    JMessage *msg = mapMsg.front();
    mapMsg.pop_front();

    rowMap[msg->x].remove(msg);
    colMap[msg->y].remove(msg);

    int BLKSIZE = msg->fsize;

    ckout << "XXXXXXXXX cblas_dgemm happening" << endl;

    cblas_dgemm(CblasRowMajor,
                CblasNoTrans, CblasNoTrans,
                BLKSIZE, BLKSIZE, BLKSIZE,
                -1.0, msg->first,
                BLKSIZE, msg->second, BLKSIZE,
                1.0, msg->LU, BLKSIZE);

    msg->block.matrixUpdated(msg->step);

    delete msg;
  }
  
  if (mapMsg.size() > 0) {
    CkEntryOptions opts;
    opts.setPriority(100000);

    thisProxy[CkMyPe()].tryAgain(0, &opts);
  }
}
    
void Scheduler::cpuFree(int cpu) {
  cpuStatus[cpu] = false;
}

void Scheduler::finishedGPU(list<JMessage> msgs) {
  for (list<JMessage>::iterator iter = msgs.begin(); 
       iter != msgs.end(); ++iter) {

    // TODO: set internal step
    iter->block.matrixUpdated(0);
  }

  GPUworking = false;
}

void Scheduler::haveWork(CProxyElement_LUBlk block, double** first, 
			 double** second, double **LU, int x, int y, int cx,
			 int cy, int fsize, int ssize, int step) {

  ckout << "haveWork called" << endl;

  msg_count++;

  JMessage *msg = new JMessage();
  //printf("*first = %d\n", **first);
  msg->first = *first;
  msg->second = *second;
  msg->LU = *LU;
  msg->block = block;
  msg->fsize = fsize;
  msg->ssize = ssize;
  msg->x = x;
  msg->y = y;
  msg->cx = cx;
  msg->cy = cy;
  msg->step = step;

  mapMsg.push_back(msg);

  rowMap[msg->x].push_back(msg);

  if (rowMap[msg->x].size() > maxRowSize) {
    maxRow = msg->x;
    maxRowSize = rowMap[msg->x].size();
  }

  colMap[msg->y].push_back(msg);

  if (colMap[msg->y].size() > maxColSize) {
    maxCol = msg->y;
    maxColSize = colMap[msg->y].size();
  }

  if (mapMsg.size() != 0) {
    CkEntryOptions opts;
    opts.setPriority(1000000);

    thisProxy[CkMyPe()].tryAgain(0, &opts);
  }
    
}

void Scheduler::checkIn() {
  //TODO: refactor so ends with this "complete" message
  //mainProxy.complete(total_size, sch_count, msg_count, cpu_count, gpu_count);
}

void Scheduler::findNewMax() {

}

list<JMessage*>* Scheduler::findLargestAgglom() {
  int elm1 = findLargestInMap(rowMap);
  int elm2 = findLargestInMap(colMap);

  int sz1 = -1, sz2 = -1;

  if (elm1 != -1)
    sz1 = rowMap[elm1].size();

  if (elm2 != -1)
    sz2 = colMap[elm2].size();

  if (sz1 == -1 && sz2 == -1)
    return NULL;
  else if (sz1 > sz2)
    return &rowMap[elm1];
  else
    return &colMap[elm2];
}

int Scheduler::findLargestInMap(map<int, list<JMessage*> > map1) {
  int maxElm = -1, maxElmSize = 0;
  bool inFirst = true;

  for (map<int, list<JMessage*> >::iterator iter = map1.begin();
       iter != map1.end(); ++iter) {
    if (iter->second.size() > maxElmSize) {
      maxElm = iter->first;
      maxElmSize = iter->second.size();
    }
  }
  
  return maxElm;
}

#include "scheduler.def.h"
