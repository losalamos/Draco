###############################################################################
# draco/src/3T/makefile
# Geoffrey Furnish
# 4 December 1997
###############################################################################
# @> Rules for constructing the 3T support library for Draco (lib3T.a).
###############################################################################

pkg := 3T

need_mesh	:= yes
need_linalg	:= yes

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

# Put non template code here

libsrcs = 

# Get the library built, minding our various other targets too.

pretarg = deps

resolve_libs = $(dep_libmesh) $(dep_liblinalg) \
	$(dep_libc4) $(dep_libnml) $(dep_libds++)

include ../config/make.lib

deps:
ifneq ($(dodeps),no)
	(cd ../ds++;   $(MAKE))
	(cd ../nml;    $(MAKE) dodeps=no)
	(cd ../c4;     $(MAKE) dodeps=no)
	(cd ../geom;   $(MAKE) dodeps=no)
	(cd ../mesh;   $(MAKE) dodeps=no)
	(cd ../linalg; $(MAKE) dodeps=no)
endif

# Rules for building docs

docsrcs	= $(wildcard *.hh *.cc)

include ../config/make.doc

# Go pick up the standard signature for makefiles.

cc_files := $(libsrcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/3T/makefile
###############################################################################
