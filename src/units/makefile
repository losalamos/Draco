###############################################################################
# draco/src/units/makefile
# Randy M. Roberts
# 7 April 1998
###############################################################################
# @> Rules for constructing the units support library for Draco (libunits.a).
###############################################################################

pkg := units

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

# Put non template code here

libsrcs = Units.cc

# Get the library built, minding our various other targets too.

pretarg = deps

resolve_libs = -L../ds++ -lds++_$(env)

include ../config/make.lib

deps:
ifneq ($(dodeps),no)
	(cd ../ds++;   $(MAKE) $(subnj))
endif

# Rules for building docs

docsrcs	= $(wildcard *.hh *.cc)

include ../config/make.doc

# Go pick up the standard signature for makefiles.

cc_files := $(libsrcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/units/makefile
###############################################################################
