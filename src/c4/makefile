###############################################################################
# gts/src/c4/makefile
# Geoffrey Furnish
# 12 January 1995
###############################################################################
# @> Rules for building the C4 library and test cases.
###############################################################################

pkg := c4

need_libm	:= yes

# Go get the standard includeable makefile for GTS components.

include ../config/makefile

# The core set of library files.

libsrcs = \
	NodeInfo.cc \
        SpinLock.cc \
	Sync.cc \
        global_nt.cc \
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

# Build the tests.

tstsrcs	:= $(wildcard t*.cc) c4run.cc

tstsrcs := $(tstsrcs:tstSort.cc=)

include ../config/make.tests

# Go pick up the standard signature for makefiles.

aux_clean := testclean

testclean:
	$(RM) t*_$(env)

cc_files	:= $(libsrcs) $(tsrcs)

include ../config/make.end

###############################################################################
# 			End of gts/src/c4/makefile
###############################################################################
