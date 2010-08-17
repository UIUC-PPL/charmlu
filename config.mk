#------------------------------------------------
# To compile on order or faith(Ubuntu):
BLAS_INC  = -DUSE_CBLAS_H=1
BLAS_LD   =
BLAS_LIBS = -lcblas

# To compile on hope (OSX):
#BLAS_INC  = -DUSE_ACCELERATE_BLAS=1
#BLAS_LD   =
#BLAS_LIBS = -framework Accelerate

# To compile against a custom atlas install on linux:
#BLAS_INC  = -I/scratch/idooley2/atlas-install/include  -DUSE_CBLAS_H=1
#BLAS_LD   = -L/scratch/idooley2/atlas-install/lib
#BLAS_LIBS = -llapack -lf77blas -lcblas -latlas


# To compile with mkl on abe
#BLAS_INC  = -DUSE_MKL_CBLAS_H=1 -I${MKL_HOME}/include/
#BLAS_LD   = -L${MKL_HOME}/lib/em64t
#BLAS_LIBS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

# To compile with ATLAS on abe
#BLAS_INC  = -I/u/ac/idooley2/LU/atlas-install/include -DUSE_CBLAS_H=1
#BLAS_LD   = -L/u/ac/idooley2/LU/atlas-install/lib
#BLAS_LIBS = -llapack -lf77blas -lcblas -latlas


# To compile on Cray XT5 with "module load atlas/3.8.3"
#BLAS_INC  = -I$(ATLAS_DIR)/include -DUSE_BLAS  -DUSE_CBLAS_H=1 -DUSE_MEMALIGN
#BLAS_LD   = -L$(ATLAS_DIR)/lib
#BLAS_LIBS = -llapack -lf77blas -lcblas -latlas

# To compile on Cray XT5 with "module load acml"
#BLAS_INC  = -I/lustre/scratch/idooley2/LU/CBLAS/include  -DUSE_CBLAS_H=1
#BLAS_LD   = -L/lustre/scratch/idooley2/LU/CBLAS/lib -L$(ACML_DIR)/gfortran64/lib/
#BLAS_LIBS = -lcblas -llibacml 


# To compile on BG/P with ESSL:
#BGP_ESSL  = /soft/apps/ESSL-4.4.1-0
#BGP_LIBS  = -L$(BGP_ESSL)/lib \
#	        -L/bgsys/ibm_compilers/sles10/prod/opt/ibmcmp/xlf/bg/11.1/bglib/ \
#	        -L/soft/apps/ibmcmp/xlsmp/bg/1.7/bglib \
#	        -L/soft/apps/ibmcmp/xlf/bg/11.1/bglib \

#BLAS_INC  = -DUSE_ESSL=1 -I$(BGP_ESSL)/include -DUSE_MEMALIGN=1
#BLAS_LD   = $(BGP_LIBS)
#BLAS_LIBS = -lesslbg -lmass -lxlfmath -lxlf90_r -lxlsmp -lxlomp_ser -lpthread


# BluePrint
#BLAS_INC = -DUSE_ESSL=1 -DUSE_MEMALIGN=1
#BLAS_LIBS = -lessl
