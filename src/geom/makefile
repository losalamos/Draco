###############################################################################
# draco/src/geom/makefile
# Geoffrey Furnish
# 31 July 1997
###############################################################################
# @> Rules for constructing the geometry library for Draco (libgeom.a).
###############################################################################

pkg := geom

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

# Put non template code here

libsrcs = \
	Cartesian1d.cc

# Get the library built, minding our various other targets too.

pretarg = deps
resolve_libs = \
	-L../c4 -lc4_$(env) \
	-L../nml -lnml_$(env) \
	-L../ds++ -lds++_$(env)

include ../config/make.lib

deps:
ifneq ($(dodeps),no)
	(cd ../ds++; $(MAKE) $(subnj))
	(cd ../nml;  $(MAKE) $(subnj) dodeps=no)
	(cd ../c4;   $(MAKE) $(subnj) dodeps=no)
endif

# Rules for building docs

docsrcs	= $(wildcard *.hh *.cc)

include ../config/make.doc

# Go pick up the standard signature for makefiles.

cc_files := $(libsrcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/geom/makefile
###############################################################################
