###############################################################################
# draco/src/radphys/makefile
# Geoffrey Furnish
# 10 July 1998
###############################################################################
# @> Rules for constructing the radphys support library for Draco (libradphys.a).
###############################################################################

pkg := radphys

need_units	:= yes

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

# Put non template code here

libsrcs = RadiationPhysics.cc

# Get the library built, minding our various other targets too.

pretarg = deps

resolve_libs = -L../units -lunits_$(env) -L../ds++ -lds++_$(env)

include ../config/make.lib

deps:
ifneq ($(dodeps),no)
	(cd ../ds++;    $(MAKE) $(subnj))
	(cd ../units;   $(MAKE) $(subnj) dodeps=no)
endif

# Rules for building docs

docsrcs	= $(wildcard *.hh *.cc)

include ../config/make.doc

# Go pick up the standard signature for makefiles.

testclean:
	(cd test; $(MAKE) clean)

aux_clean := testclean

cc_files := $(libsrcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/radphys/makefile
###############################################################################
