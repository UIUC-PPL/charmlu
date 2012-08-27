include config.mk
#include ../config.mk

# The relevant source files for this project
RAWSRC    = lu.C scheduler.C benchmark.C
INTF      = luUtils.ci lu.ci driver.ci

# Specify the exe name and the arguments to run it with
NP        = 4
TARGET    = charmlu
BINS      = $(TARGET)
ARGS      = 64 16 500 8 2

# Specify the compilers, run script, flags etc.
OPT       = -O0
#CPPFLAGS += -DSCHED_PIVOT_REDN
#CPPFLAGS += -DCHARMLU_USEG_FROM_BELOW
CPPFLAGS += -DSEND_LIM=$(SEND_LIM) $(BLAS_INC)
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

REVNUM = 
ifneq ($(REVNUM),)
  CPPFLAGS += -DLU_REVISION=$(REVNUM)
endif

.PHONY: all clean realclean again test bgtest translateInterface

all: $(BINS)

charmlu: CXX = $(CHARMC)
charmlu: $(RAWSRC:.C=.o)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	$(RM) $(wildcard *.decl.h *.def.h *.d *.di *.o *.stamp) charmrun

realclean: clean
	$(RM) $(wildcard *.log *.projrc *.sts) $(BINS)

again: 
	$(MAKE) clean; $(MAKE)

test: all
	@echo "########################################################################################"
	$(EXEC) $(EXECFLAGS) $(TARGET) $(ARGS)

regtest: $(TARGET)
	perl $(TEST_SCRIPT) $(EXEC) $(EXECFLAGS) $(TARGET)

# A test program for getting essl working on BG/P
#link-test-essl : link-test-essl.cxx
#	/soft/apps/ibmcmp-jan2010/vacpp/bg/9.0/bin/bgxlc++ link-test-essl.cxx $(BLAS_INC) $(BLAS_LD)

CHARMINC = $(shell dirname $(CHARMC))/../include
CHARMBIN = $(shell dirname $(CHARMC))

####### Pattern rules
# Rule to generate dependency information for C++ source files
%.d: %.C
	$(info Generating dependencies for $<)
	@g++ -MM -MG $(CPPFLAGS) $(INCDIRS:%=-I%) -I$(CHARMINC) $< | perl $(CHARMBIN)/dep.pl $(CHARMINC) > $@
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
-include $(RAWSRC:.C=.d)
-include $(INTF:.ci=.di)
endif
endif

