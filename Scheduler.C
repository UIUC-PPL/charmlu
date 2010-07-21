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

  //ckout << "calling tryAgain" << endl;

  // GPU offloading
  if (1) {
    int totalSize = 0;
    int numberAgglom = 0;
    list<JMessage*> toOffload;
    list<JMessage*> toRemove;

    if (mapMsg.size() >= MIN_SIZE) {
      while (numberAgglom < MAX_SIZE) {
        list<JMessage*>* msgs = findLargestAgglom();
        bool cont = true;

        //ckout << "doing this" << endl;

        if (msgs != NULL) {
          for (list<JMessage*>::iterator iter = msgs->begin();
               iter != msgs->end() && numberAgglom < MAX_SIZE; ++iter) {
            int tsize = (*iter)->fsize + (*iter)->ssize;

            if (tsize + totalSize < MAX_OFFLOAD_SIZE) {
              toOffload.push_back(*iter);
              toRemove.push_back(*iter);

              numberAgglom++;
            } else {
              cont = false;
              break;
            }
          }
        } else {
          cont = false;
          //break;
        }

        for (list<JMessage*>::iterator iter = toRemove.begin();
             iter != toRemove.end(); ++iter) {
          removeItemXY(mapMsg, *iter);
          removeItemXY(rowMap[(*iter)->x], *iter);
          removeItemXY(colMap[(*iter)->y], *iter);
        }

        toRemove.clear();

        if (!cont)
          break;
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
  if (0 && mapMsg.size() > 0) {
    cpu_count++;

    JMessage *msg = mapMsg.front();
    mapMsg.pop_front();

    rowMap[msg->x].remove(msg);
    colMap[msg->y].remove(msg);

    int BLKSIZE = msg->fsize;

    square_dgemm(BLKSIZE, msg->second, msg->first, msg->LU);

    for (int i = 0; i < msg->fsize*msg->ssize; i++) {
      printf("CPU work: msg->LU[%d] = %f\n", i, msg->LU[i]);
    }

    /*cblas_dgemm(CblasRowMajor,
                CblasNoTrans, CblasNoTrans,
                BLKSIZE, BLKSIZE, BLKSIZE,
                -1.0, msg->first,
                BLKSIZE, msg->second, BLKSIZE,
                1.0, msg->LU, BLKSIZE);*/

    msg->block.matrixUpdated(msg->step);

    delete msg;
  }
  
  if (mapMsg.size() > 0) {
    CkEntryOptions opts;
    opts.setPriority(100000);

    thisProxy[CkMyPe()].tryAgain(0, &opts);
  }
}

void Scheduler::square_dgemm (int M, double *A, double *B, double *C) {
  int i, j, k;

  for (i = 0; i < M*M; i++) {
    A[i] = A[i] * -1.0;
  }

  for (i = 0; i < M; ++i) {
    for (j = 0; j < M; ++j) {
      const double *Ai_ = A + i;
      const double *B_j = B + j*M;

      double cij = *(C + j*M + i);

      for (k = 0; k < M; ++k) {
        cij += *(Ai_ + k*M) * *(B_j + k);
      }

      *(C + j*M + i) = cij;
    }
  }
}

void Scheduler::cpuFree(int cpu) {
  cpuStatus[cpu] = false;
}

void Scheduler::finishedGPU(list<JMessage> msgs) {
  for (list<JMessage>::iterator iter = msgs.begin(); 
       iter != msgs.end(); ++iter) {

    //square_dgemm(iter->fsize, iter->second, iter->first, iter->LU);

    iter->block.matrixUpdated(iter->step);
  }

  GPUworking = false;
}

void Scheduler::haveWork(CProxyElement_LUBlk block, double** first, 
			 double** second, double **LU, int x, int y, int cx,
			 int cy, int fsize, int ssize, int step) {

  //ckout << "haveWork called" << endl;

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

  /*for (int i = 0; i < fsize * ssize; i++) {
    ckout << "L = " << msg->first[i] << endl;
    ckout << "U = " << msg->second[i] << endl;
    ckout << "LU = " << msg->LU[i] << endl;
    }*/

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
  mainProxy.complete(total_size, sch_count, msg_count, cpu_count, gpu_count);
}

list<JMessage*>* Scheduler::findLargestAgglom() {
  int elm1 = findLargestInMap(rowMap);
  int elm2 = findLargestInMap(colMap);

  //ckout << "--------------elm1 = " << elm1 << ", " << "elm2 = " << elm2 << endl;

  int sz1 = -1, sz2 = -1;

  if (elm1 != -1)
    sz1 = rowMap[elm1].size();

  if (elm2 != -1)
    sz2 = colMap[elm2].size();

  //ckout << "--------------sz1 = " << sz1 << ", " << "sz2 = " << sz2 << endl;

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

void Scheduler::removeItemXY(list<JMessage*>& list1, JMessage* &msg) {
  for (list<JMessage*>::iterator iter = list1.begin();
       iter != list1.end(); ++iter) {
    if ((*iter)->x == msg->x &&
        (*iter)->y == msg->y &&
        (*iter)->step == msg->step) {
      list1.erase(iter);
    }
  }
}

#include "scheduler.def.h"
