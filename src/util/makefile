###############################################################################
# draco/src/util/makefile
# Geoffrey Furnish
# 24 September 1997
###############################################################################
# @> Rules for constructing the utilities library for Draco.
###############################################################################

pkg := util

need_util	:= yes

# Go get the standard includeable makefile for GTS components.

include ../config/makefile

libsrcs = \
	ADKey.cc \
	ADFile.cc \
	dump_nt.cc \
	FileSet.cc \
	Timer.cc \
	persist.cc \
	random.cc

pretarg = deps
resolve_libs = \
	-L../spdf -lspdf_$(env) \
	-L../c4 -lc4_$(env) \
	-L../ds++ -lds++_$(env)

include ../config/make.lib

deps:
ifneq ($(dodeps),no)
	(cd ../ds++; $(MAKE) $(subnj))
	(cd ../spdf; $(MAKE) $(subnj) dodeps=no)
endif

# Now concern ourselves with the business of test programs.

notstlib := 1
tstsrcs	:= $(wildcard tst*.cc)

include ../config/make.tests

# Go pick up the standard signature for makefiles.

cc_files := $(libsrcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/util/makefile
###############################################################################
