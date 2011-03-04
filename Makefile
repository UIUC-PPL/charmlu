include config.mk

# Point this to your charm installation
CHARM ?= $(HOME)/charm

# Charm directory structure
CHARMBIN := $(CHARM)/bin
CHARMINC := $(CHARM)/include

# The relevant source files for this project
HDR       =
SRC       = lu.C scheduler.C
OBJ       = $(SRC:.C=.o)
INTF      = lu.ci   

# Specify the exe name and the arguments to run it with
NP        = 4
TARGET    = lu
ARGS      = 64 16 500 8 2

# Specify the compilers, run script, flags etc.
CXX       = $(CHARMBIN)/charmc
CHARMC    = $(CHARMBIN)/charmc
OPT       = -g -O3
CPPFLAGS += -DADAPT_SCHED_MEM $(BLAS_INC)
CXXFLAGS += -language charm++ $(OPT)
LDFLAGS  += -module comlib -module CkMulticast -tracemode projections $(BLAS_LD)
LDLIBS   += $(BLAS_LIBS)
EXEC      = ./charmrun
EXECFLAGS = +p$(NP) ++local
TEST_SCRIPT=test_lu.pl
ifdef $(NODELIST)
  EXECFLAGS += ++nodelist $(NODELIST)
endif


########### This stuff should be able take care of itself ############

.PHONY: all clean realclean again test bgtest translateInterface

all: $(TARGET)

$(TARGET): $(OBJ) 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	$(RM) $(wildcard *.decl.h *.def.h *.d *.di *.o *.stamp) $(TARGET) charmrun

realclean: clean
	$(RM) $(wildcard *.log *.projrc *.sts)

again: 
	$(MAKE) clean; $(MAKE)

test: all
	@echo "########################################################################################"
	$(EXEC) $(EXECFLAGS) $(TARGET) $(ARGS)

regtest: all
	perl $(TEST_SCRIPT) $(EXEC) $(EXECFLAGS) $(TARGET)

# A test program for getting essl working on BG/P
#link-test-essl : link-test-essl.cxx
#	/soft/apps/ibmcmp-jan2010/vacpp/bg/9.0/bin/bgxlc++ link-test-essl.cxx $(BLAS_INC) $(BLAS_LD)

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
	$q$(CHARMC) $< && touch $@

# Include the generated files containing dependency info
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),realclean)
-include $(SRC:.C=.d)
-include $(INTF:.ci=.di)
endif
endif

