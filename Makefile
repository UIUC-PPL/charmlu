OPTS=-g -O0
#OPTS = -O3



#------------------------------------------------
# To compile on order or faith(Ubuntu):
#BLAS_INC = -DUSE_BLAS -DUSE_CBLAS_H
#BLAS_LD = -lcblas

# To compile on hope (OSX):
#BLAS_INC = -DUSE_BLAS -DUSE_ACCELERATE_BLAS=1
#BLAS_LD = -framework Accelerate

# To compile against a custom atlas install on linux:
BLAS_INC = -I/scratch/idooley2/atlas-install -DUSE_BLAS  -DUSE_CBLAS_H
BLAS_LD = -L/scratch/idooley2/atlas-install -lcblas -latlas


# ----------------------------------------------
# Some options that have been used in the past

#BLAS_LD =  -lm /expand/home/idooley2/src/ATLAS/Linux_insight/lib/libatlas.a /expand/home/idooley2/src/ATLAS/Linux_insight/lib/libcblas.a
#BLAS_INC = -I/expand/home/idooley2/local/blas/include -DUSE_BLAS  
#BLAS_LD  =  -L/expand/home/idooley2/local/lib -lm -lcblas -latlas
#BLAS_LD = -L$(MKL_HOME)/lib/em64t -lmkl_lapack -lmkl_em64t -lmkl_core -lpthread
#BLAS_LD = $(MKL_HOME)/lib/em64t/libmkl_intel_lp64.a $(MKL_HOME)/lib/em64t/libmkl_sequential.a $(MKL_HOME)/lib/em64t/libmkl_core.a -lguide -lpthread -lm
#BLAS_INC = -DUSE_BLAS -I$(MKL_HOME)/include





PROJ = -tracemode summary -tracemode projections
#MULTICAST = -module CkMulticast

CHARMC=../charm/bin/charmc $(OPTS) 
#CHARMC=${HOME}/current/charm/net-linux/bin/charmc $(OPTS) -g
#CHARMC=${HOME}/current/lastestfromcvs/charm/net-linux-amd64/bin/charmc $(OPTS)
#CHARMC=${HOME}/charm/bin/charmc $(FLAGS)
#CHARMC=${HOME}/current/new/charm/mpi-linux-amd64-vmi-mpicxx/bin/charmc $(OPTS)
#CHARMC=${HOME}/current/new/charm/mpi-linux-amd64-smp-vmi-mpicxx/bin/charmc $(OPTS)
#CHARMC=${HOME}/current/new/charm/net-linux-amd64-smp-icc10/bin/charmc $(OPTS)
#CHARMC=${HOME}/current/new/charm/net-linux-amd64-icc10/bin/charmc $(OPTS)
#CHARMC=charm/bin/charmc $(OPTS)


MODULES=  -module ControlPoints  -tracemode controlPoints -tracemode projections

#MEM= -memory gnu


all: lu lu-proj

lu: lu.o 
	$(CHARMC) -language charm++ -o lu lu.o $(BLAS_LD) $(MULTICAST)  $(MODULES) $(MEM)

lu-proj: lu.o 
	$(CHARMC) -language charm++ -o lu-proj lu.o $(BLAS_LD) $(PROJ) $(MULTICAST)  $(MODULES) $(MEM)

lu.decl.h: lu.ci
	$(CHARMC)  lu.ci

clean:
	rm -f *.decl.h *.def.h conv-host *.o charmrun *~ lu lu-blas lu-mem lu-blas-proj.*.log lu-blas-proj.*.sum lu-blas-proj.*.sts lu-blas-proj.sts lu-blas-proj.projrc lu-blas-proj lu-proj controlPointData.txt lu*.log lu*.sum lu*.sts lu*.projrc SummaryDump.out

lu.o: lu.C lu.decl.h
	$(CHARMC) -c lu.C -o lu.o $(BLAS_INC) $(OPTS)
