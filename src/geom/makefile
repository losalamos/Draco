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
	../c4/libc4_$(env).$(lib_suffix) \
	../nml/libnml_$(env).$(lib_suffix) \
	../ds++/libds++_$(env).$(lib_suffix)

include ../config/make.lib

deps:
ifneq ($(dodeps),no)
	(cd ../ds++; $(MAKE))
	(cd ../nml;  $(MAKE) dodeps=no)
	(cd ../c4;   $(MAKE) dodeps=no)
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
