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

resolve_libs = \
	-L../c4 -lc4_$(env) \
	-L../nml -lnml_$(env) \
	-L../ds++ -lds++_$(env)

include ../config/make.lib

# Test programs

tstsrcs := $(wildcard t*.cc)

tstpretarg := deps

include ../config/make.tests

# Get DejaGnu support

include ../config/make.deja

# Dependencies

deps:
ifneq ($(dodeps),no)
	(cd ../ds++;   $(MAKE) $(subnj))
	(cd ../nml;    $(MAKE) $(subnj) dodeps=no)
	(cd ../c4;     $(MAKE) $(subnj) dodeps=no)
	(cd ../geom;   $(MAKE) $(subnj) dodeps=no)
endif

# Rules for building docs

docsrcs	= $(wildcard *.hh *.cc)

include ../config/make.doc

# Go pick up the standard signature for makefiles.

cc_files := $(libsrcs) $(tstsrcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/mesh/makefile
###############################################################################
