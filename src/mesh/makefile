###############################################################################
# draco/src/mesh/makefile
# Geoffrey Furnish
# 13 June 1997
###############################################################################
# @> Rules for constructing the mesh library for Draco (libmesh.a).
###############################################################################

pkg := mesh

need_mesh	:= yes

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

# Put non template code here

libsrcs = \
	Mesh_DB.cc Mesh_XYZ.cc

#	K_Mesh_DB.cc Mesh_DB.cc

# Get the library built, minding our various other targets too.

pretarg = deps

resolve_libs = $(dep_libc4) $(dep_libnml) $(dep_libds++) #$(dep_libpooma)


include ../config/make.lib

# Test programs

tstsrcs := $(wildcard t*.cc)

tstpretarg := deps

include ../config/make.tests

# Testing out some dejagnu support

deja:
	(cd testsuites; $(MAKE) links)
	runtest --tool=ds++ --srcdir=./testsuites ENV=$(env) DBC=$(DBC)

# Dependencies

deps:
ifneq ($(dodeps),no)
	(cd ../ds++;   $(MAKE))
	(cd ../nml;    $(MAKE) dodeps=no)
	(cd ../c4;     $(MAKE) dodeps=no)
	(cd ../geom;   $(MAKE) dodeps=no)
endif

# Rules for building docs

docsrcs	= $(wildcard *.hh *.cc)

include ../config/make.doc

# Go pick up the standard signature for makefiles.

cc_files := $(libsrcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/mesh/makefile
###############################################################################
