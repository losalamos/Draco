###############################################################################
# gts/src/spdf/makefile
# Geoffrey Furnish
# 11 August 1995
###############################################################################
# @> Rules for building the SPDF library and test cases.
###############################################################################

pkg := spdf

# Go get the standard includeable makefile for GTS components.

include ../config/makefile

# The core set of library files.

libsrcs = 

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

docsrcs	= $(wildcard *.h *.cc)

include ../config/make.doc

# Rules for building test executables

tstsrcs		= $(wildcard t*.cc)

include ../config/make.tests

# Go pick up the standard signature for makefiles.

cc_files := $(libsrcs) $(tstsrcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/spdf/makefile
###############################################################################
