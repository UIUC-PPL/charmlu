/*
 * Scheduler.C
 *
 *  Created on: July 14, 2010
 *      Author: Jonathan
 */

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

  if (1 || !GPUworking) {
    int totalParticles = 0;
    int numberAgglom = 0;
    list<JMessage*> largest;
    list<JMessage*>::iterator iter;

    for (iter = mapMsg.begin(); iter != mapMsg.end(); ++iter) {
      if (totalParticles + (*iter)->getSize() > MAX_OFFLOAD_SIZE
	  || numberAgglom == MAX_SIZE)
	break;

      totalParticles += (*iter)->getSize();
      numberAgglom++;
    }

    if (numberAgglom >= MIN_SIZE) {
      //ckout << "before splice " << mapMsg.size();

      //largest.splice(largest.begin(), mapMsg, mapMsg.begin(), iter);

      for (list<JMessage*>::iterator iter2 = mapMsg.begin(); 
	   iter2 != iter; iter2++) {
	largest.push_back(*iter2);
      }

      mapMsg.erase(mapMsg.begin(), iter);

      int size = largest.size();

      //ckout << "numberAgglom = " << numberAgglom << ", queue size = " << mapMsg.size() << endl; 

      //ckout << ", splice happened mapMsg = " << mapMsg.size() << ", largest = " << largest.size() << endl;

      // run on GPU
      gpu_count++;
      total_size += size;

      int partsize = 0, k = 0;
	
      list<JMessage> toSend;

      for (list<JMessage*>::iterator iter2 = largest.begin(); 
	   iter2 != largest.end(); ++iter2) {
	toSend.push_back(**iter2);
      }

      GPUworking = true;

      CkEntryOptions opts;
      opts.setPriority(-1000);

      //ckout << "running " << numberAgglom  << " msg on GPU, " << "queue size is " << mapMsg.size() << endl;

      //gpu_state++;
      //gpu_state = gpu_state % 4;

      gpu[CkMyPe()].gpu_offload(toSend, &opts);
    }
  }

  //for (int = 0; i < cpuStatus.size(); i++) {
  if (1 || !cpuStatus[thisIndex]) {
    //ckout << "Found free cpu" << endl;

    cpuStatus[thisIndex] = true;

    //ckout << "running " << (MAX_CPU_SIZE > mapMsg.size() ? mapMsg.size() : MAX_CPU_SIZE) << " msg on CPU, " << "queue size is " << mapMsg.size() << endl;
	
    // there is work to do: run one message on CPU
    for (int i = 0; i < MAX_CPU_SIZE && mapMsg.size() > 0; i++) {
      cpu_count++;

      JMessage *msg = mapMsg.front();
      mapMsg.pop_front();

      //TODO: fix this CPU work
      //msg->block.matrixUpdated(0);

      delete msg;
    }
  }
  //}

  if (mapMsg.size() > 0) {
    CkEntryOptions opts;
    opts.setPriority(100000);

    thisProxy[thisIndex].tryAgain(0, &opts);
  }
}
    
void Scheduler::cpuFree(int cpu) {
  /*char note[200];

    sprintf(note, "cpu free happening");
    traceUserSuppliedNote(note);*/

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
			 double** second, int x, int y, int cx,
			 int cy, int fsize, int ssize) {

  //ckout << "haveWork called" << endl;

  msg_count++;

  JMessage *msg = new JMessage();
  //printf("*first = %d\n", **first);
  msg->first = *first;
  msg->second = *second;
  msg->block = block;
  msg->fsize = fsize;
  msg->ssize = ssize;
  msg->x = x;
  msg->y = y;
  msg->cx = cx;
  msg->cy = cy;

  mapMsg.push_back(msg);

  if (mapMsg.size() != 0) {
    CkEntryOptions opts;
    opts.setPriority(1000000);

    thisProxy[thisIndex].tryAgain(0, &opts);
  }
    
}

void Scheduler::checkIn() {
  //TODO: refactor so ends with this "complete" message
  //mainProxy.complete(total_size, sch_count, msg_count, cpu_count, gpu_count);
}

#include "scheduler.def.h"
