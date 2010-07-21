##OPTS=-O3 -Q -qhot -qarch=450d -qtune=450
OPTS=-g -O0

#------------------------------------------------
# To compile on order or faith(Ubuntu):
#BLAS_INC =  -DUSE_CBLAS_H=1
#BLAS_LD = -lcblas

# To compile on hope (OSX):
BLAS_INC = -DUSE_ACCELERATE_BLAS=1
BLAS_LD = -framework Accelerate

# To compile against a custom atlas install on linux:
#BLAS_INC = -I/scratch/idooley2/atlas-install/include  -DUSE_CBLAS_H=1
#BLAS_LD = -L/scratch/idooley2/atlas-install/lib -llapack -lf77blas -lcblas -latlas


# To compile with mkl on abe
#BLAS_INC = -DUSE_MKL_CBLAS_H=1 -I${MKL_HOME}/include/
#BLAS_LD = -L${MKL_HOME}/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

# To compile with ATLAS on abe
#BLAS_INC = -I/u/ac/idooley2/LU/atlas-install/include -DUSE_CBLAS_H=1
#BLAS_LD = -L/u/ac/idooley2/LU/atlas-install/lib -llapack -lf77blas -lcblas -latlas


# To compile on Cray XT5 with "module load atlas/3.8.3"
#BLAS_LD = -L$(ATLAS_DIR)/lib -llapack -lf77blas -lcblas -latlas
#BLAS_INC = -I$(ATLAS_DIR)/include -DUSE_BLAS  -DUSE_CBLAS_H=1 -DUSE_MEMALIGN

# To compile on Cray XT5 with "module load acml"
#BLAS_LD = -L/lustre/scratch/idooley2/LU/CBLAS/lib -lcblas -L$(ACML_DIR)/gfortran64/lib/ -llibacml 
#BLAS_INC = -I/lustre/scratch/idooley2/LU/CBLAS/include  -DUSE_CBLAS_H=1


# To compile on BG/P with ESSL:
#BGP_ESSL = /soft/apps/ESSL-4.4.1-0
#BLAS_INC = -DUSE_ESSL=1 -I$(BGP_ESSL)/include -DUSE_MEMALIGN=1
#BGP_LIBS = -L$(BGP_ESSL)/lib \
#	-L/bgsys/ibm_compilers/sles10/prod/opt/ibmcmp/xlf/bg/11.1/bglib/ \
#	-L/soft/apps/ibmcmp/xlsmp/bg/1.7/bglib \
#	-L/soft/apps/ibmcmp/xlf/bg/11.1/bglib \
#	 -lesslbg -lmass   -lxlfmath -lxlf90_r -lxlsmp -lxlomp_ser -lpthread

#BLAS_LD =  $(BGP_LIBS)


# A test program for getting essl working on BG/P
#link-test-essl : link-test-essl.cxx
#	/soft/apps/ibmcmp-jan2010/vacpp/bg/9.0/bin/bgxlc++ link-test-essl.cxx $(BLAS_INC) $(BLAS_LD)



# ----------------------------------------------
# Some options that have been used in the past

#BLAS_LD =  -lm /expand/home/idooley2/src/ATLAS/Linux_insight/lib/libatlas.a /expand/home/idooley2/src/ATLAS/Linux_insight/lib/libcblas.a
#BLAS_INC = -I/expand/home/idooley2/local/blas/include -DUSE_BLAS  
#BLAS_LD  =  -L/expand/home/idooley2/local/lib -lm -lcblas -latlas
#BLAS_LD = -L$(MKL_HOME)/lib/em64t -lmkl_lapack -lmkl_em64t -lmkl_core -lpthread
#BLAS_LD = $(MKL_HOME)/lib/em64t/libmkl_intel_lp64.a $(MKL_HOME)/lib/em64t/libmkl_sequential.a $(MKL_HOME)/lib/em64t/libmkl_core.a -lguide -lpthread -lm
#BLAS_INC = -DUSE_BLAS -I$(MKL_HOME)/include



PROJ = -tracemode projections
#MULTICAST = -module CkMulticast

CHARMC=charmc $(OPTS)
#CHARMC=${HOME}/current/charm/net-linux/bin/charmc $(OPTS) -g
#CHARMC=${HOME}/current/lastestfromcvs/charm/net-linux-amd64/bin/charmc $(OPTS)
#CHARMC=${HOME}/charm/bin/charmc $(FLAGS)
#CHARMC=${HOME}/current/new/charm/mpi-linux-amd64-vmi-mpicxx/bin/charmc $(OPTS)
#CHARMC=${HOME}/current/new/charm/mpi-linux-amd64-smp-vmi-mpicxx/bin/charmc $(OPTS)
#CHARMC=${HOME}/current/new/charm/net-linux-amd64-smp-icc10/bin/charmc $(OPTS)
#CHARMC=${HOME}/current/new/charm/net-linux-amd64-icc10/bin/charmc $(OPTS)
#CHARMC=charm/bin/charmc $(OPTS)

MODULES=-module comlib

all: lu-proj

lu: lu.o Scheduler.o gpu_offloader.o gpu_kernel.o
	$(CHARMC) -language charm++ -o $@ $^  $(BLAS_LD) $(MULTICAST)  $(MODULES)

lu-proj: lu.o Scheduler.o gpu_offloader.o gpu_kernel.o
	$(CHARMC) -L/usr/local/cuda/lib -lcublas -lcudart -language charm++ -o $@ $^ $(BLAS_LD) $(MULTICAST) $(MODULES) -DADAPT_SCHED_MEM

lu.decl.h: lu.ci
	$(CHARMC)  lu.ci -DADAPT_SCHED_MEM

clean:
	rm -f *.decl.h *.def.h conv-host *.o charmrun *~ lu lu-blas lu-mem lu-blas-proj.*.log lu-blas-proj.*.sum lu-blas-proj.*.sts lu-blas-proj.sts lu-blas-proj.projrc lu-blas-proj lu-proj controlPointData.txt lu*.log lu*.sum lu*.sts lu*.projrc SummaryDump.out *.output *.error *.cobaltlog traces/* core.* perfCounterBGP.o job-lu-* moduleinit* moduleInit*

lu.o: lu.C lu.decl.h Scheduler.o gpu_offloader.o gpu_kernel.o c_common.h$
	$(CHARMC) -c lu.C -o lu.o $(BLAS_INC) $(OPTS) -DADAPT_SCHED_MEM

Scheduler.o: Scheduler.C scheduler.decl.h c_common.h$
	$(CHARMC) -c Scheduler.C -o $@ $(BLAS_INC) $(OPTS) -DADAPT_SCHED_MEM

gpu_offloader.o: gpu_offloader.C gpuwork.decl.h c_common.h$
	$(CHARMC) -c gpu_offloader.C -o $@ $(BLAS_INC) $(OPTS) -DADAPT_SCHED_MEM

NVCC = nvcc

gpu_kernel.o : gpu_kernel.cu c_common.h
	$(NVCC) $(OPTS1) -c $<

 # run for up to 15 minutes on 16 nodes * 4 pe/node. Matrix size 8192*8192
run-BGP: lu-proj
	rm -fr traces
	mkdir traces
	qsub -n 64 --mode smp -t 30 ./lu-proj 32768 1000 2 +CPSaveData +CPExhaustiveSearch +CPLoadData +traceroot traces 
#+logsize 10000000b

run: lu-proj
	-rm -rf traces
	mkdir traces
	charmrun +p2 ./lu-proj 64 1000 1 8
#++debug-no-pause
#+CPSaveData +traceroot traces +logsize 10000000