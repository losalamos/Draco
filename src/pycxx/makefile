###############################################################################
# draco/src/pycxx/makefile
# Geoffrey Furnish
# 2 April 1998
###############################################################################
# @> Rules for building the Python C++ support library
###############################################################################

pkg := pycxx

need_python	:= yes

# Go get the standard includeable makefile for components.

include ../config/makefile

# The core set of library files.

libsrcs = \
	cxxmodule.cc cxxextensions.cc

# Get the library built, minding our various other targets too.

#pretarg = deps
#auxtarg = $(TEST_PROGS)
#resolve_libs = -L../ds++ -lds++_$(env)

include ../config/make.lib

deps:
ifneq ($(dodeps),no)
	(cd ../ds++; $(MAKE) $(subnj) "do=")
endif

# Go pick up the standard signature for makefiles.

cc_files	:= $(libsrcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/pycxx/makefile
###############################################################################
