###############################################################################
# gts/src/c4/makefile
# Geoffrey Furnish
# 12 January 1995
###############################################################################
# @> Rules for building the C4 library and test cases.
###############################################################################

pkg := c4

# Go get the standard includeable makefile for GTS components.

include ../config/makefile

# The core set of library files.

libsrcs = \
	NodeInfo.cc \
        SpinLock.cc \
	Sync.cc \
        global.cc \
	C4_Req.cc

#	DeepBSwap.cc
#	DeepObj.cc \

# Pathetic totalview can't handle this file, grrrrr.
#        global.cc \
#        global_shmem.cc \

# Template definition files, for C-front systems.

tmplsrcs	= \
	BSwap.cc \
	Baton.cc \
	ArraySwap.cc \
	DstArray.cc \
	DynArraySwap.cc

# Build specified test executables using the syntax:
#	"make cxx=<env> do=<program name>".

ifdef do
TEST_PROGS = $(do:%=%_$(env))
endif

# Get the library built, minding our various other targets too.

pretarg = deps
auxtarg = $(TEST_PROGS)

include ../config/make.lib

deps:
ifneq ($(dodeps),no)
	(cd ../ds++; $(MAKE) $(subnj) "do=")
endif

# Rules for building docs

docsrcs := $(wildcard *.h *.cc)

include ../config/make.doc

# Rules for building test executables

tsrcs		= $(wildcard t*.cc)
testprogs_all	= $(tsrcs:%.cc=%_$(env))

tests: libc4_$(env).$(lib_suffix) $(testprogs_all)

$(testprogs_all) : t%_$(env) : $(env) $(env)/t%.o libc4_$(env).$(lib_suffix)
	$(CXX) $(CXXFLAGS) -o $@ $(env)/t$*.o -L. -lc4_$(env) $(LINK_LIBS) -lm

# Go pick up the standard signature for makefiles.

aux_clean := testclean

testclean:
	$(RM) t*_$(env)

cc_files	:= $(libsrcs) $(tsrcs)

include ../config/make.end

###############################################################################
# 			End of gts/src/c4/makefile
###############################################################################
