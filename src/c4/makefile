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
	C4_Req.cc

ifeq ($(threads),yes)
libsrcs += \
	ThreadGroupMember.cc \
	ThreadControl.cc
endif

# Add in the platform specific files

ifeq ($(c4),s2)
libsrcs += shmem_sr.cc shm_reduce.cc
endif
ifeq ($(c4),shmem)
libsrcs += global_shmem.cc shm_reduce.cc
endif
ifeq ($(c4),mpi)
libsrcs	+= global_mpi.cc
endif
ifeq ($(c4),no)
libsrcs	+= global_scalar.cc
endif

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
resolve_libs = -L../ds++ -lds++_$(env)

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

# Exclude shmem-rewrite specific tests if not using the shmem rewrite.
ifneq ($(c4),s2)
tstsrcs := $(tstsrcs:tsr2.cc=)
endif

# Exclude the threaded tests if not using threads for this build.
ifneq ($(threads),yes)
tstsrcs := $(tstsrcs:tthrd1.cc=)
endif

tstpretarg := $(pretarg)

include ../config/make.tests

# Testing out some dejagnu support

include ../config/make.deja

# Go pick up the standard signature for makefiles.

aux_clean := testclean

testclean:
	$(RM) t*_$(env)

cc_files	:= $(libsrcs) $(tstsrcs)

include ../config/make.end

###############################################################################
# 			End of gts/src/c4/makefile
###############################################################################
