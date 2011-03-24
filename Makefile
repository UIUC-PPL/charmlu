include config.mk

# The relevant source files for this project
RAWSRC    = lu.C scheduler.C
INTF      = lu.ci

BENCHSRC = dgerBenchmark.C
BENCHCI = dgerBenchmark.ci

# Specify the exe name and the arguments to run it with
NP        = 4
TARGET    = lu.prod
BINS      = lu.prod lu.trace
BENCH     = lu_dger
ARGS      = 64 16 500 8 2

# Specify the compilers, run script, flags etc.
CHARMC    = $(CHARMPROD)/bin/charmc
CHARMINC  = $(CHARMPROD)/include
OPT       = -O3
CPPFLAGS += -DADAPT_SCHED_MEM -DSEND_LIM=$(SEND_LIM) $(BLAS_INC)
CXXFLAGS += -language charm++ $(OPT)
LDFLAGS  += -module comlib -module CkMulticast $(BLAS_LD)
LDLIBS   += $(BLAS_LIBS)
EXEC      = ./charmrun
EXECFLAGS = +p$(NP) ++local
TEST_SCRIPT=test_lu.pl
ifdef $(NODELIST)
  EXECFLAGS += ++nodelist $(NODELIST)
endif


########### This stuff should be able take care of itself ############

# The base directory of this project
base     ?= .
GIT       = $(shell which git)
# Compute the revision number (hash) of the build and feed it to the code
ifeq ($(GIT),)
  $(warning Cannot find the git binary. Will not compute the revision number)
else
  REVNUM  = $(shell $(GIT) --git-dir=$(base)/.git rev-parse HEAD)
endif
ifneq ($(REVNUM),)
  CPPFLAGS += -DLU_REVISION=$(REVNUM)
endif


.PHONY: all clean realclean again test bgtest translateInterface

all: $(BINS)

lu_dger: CXX = $(CHARMPROD)/bin/charmc
lu_dger: $(BENCHSRC:.C=.o)
	$(CHARMC) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

lu.prod: CXX = $(CHARMPROD)/bin/charmc
lu.prod: $(RAWSRC:.C=-prod.o)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

lu.trace: CXX = $(CHARMTRACE)/bin/charmc
lu.trace: CPPFLAGS += -DLU_TRACING
lu.trace: $(RAWSRC:.C=-trace.o)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS) -tracemode projections

clean:
	$(RM) $(wildcard *.decl.h *.def.h *.d *.di *.o *.stamp) charmrun

realclean: clean
	$(RM) $(wildcard *.log *.projrc *.sts) $(BINS)

again: 
	$(MAKE) clean; $(MAKE)

test: all
	@echo "########################################################################################"
	$(EXEC) $(EXECFLAGS) $(TARGET) $(ARGS)

regtest: lu.prod
	perl $(TEST_SCRIPT) $(EXEC) $(EXECFLAGS) $(TARGET)

# A test program for getting essl working on BG/P
#link-test-essl : link-test-essl.cxx
#	/soft/apps/ibmcmp-jan2010/vacpp/bg/9.0/bin/bgxlc++ link-test-essl.cxx $(BLAS_INC) $(BLAS_LD)

####### Pattern rules
# Rule to generate dependency information for C++ source files
%.d: %.C
	$(info Generating dependencies for $<)
	@g++ -MM -MG $(CPPFLAGS) $(INCDIRS:%=-I%) -I$(CHARMINC) $< | perl $(CHARMPROD)/bin/dep.pl $(CHARMINC) > $@
#	@$(SHELL) -ec 'g++ -MM -MG $(CPPFLAGS) $(INCDIRS:%=-I%) $< \
#	| sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
#	[ -s $@ ] || rm -f $@'

# Rule to generate dependency info for charm++ interface (ci) definition files
%.di: %.ci
	$(info Generating dependencies for $<)
	@$(CHARMC) -M $< > $@

# Compilation rule for ci files
%.ci.stamp: %.ci
	$(info-ci)
	$q$(CHARMC) $(CPPFLAGS) -E $< && touch $@

# Include the generated files containing dependency info
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),realclean)
-include $(RAWSRC:.C=-prod.d) $(RAWSRC:.C=-trace.d) $(BENCHSRC:.C=.d)
-include $(INTF:.ci=.di) $(BENCHCI:.ci=.di)
endif
endif

