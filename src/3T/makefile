###############################################################################
# draco/src/3T/makefile
# Geoffrey Furnish
# 4 December 1997
###############################################################################
# @> Rules for constructing the 3T support library for Draco (lib3T.a).
###############################################################################

pkg := 3T

need_units	:= yes
need_mesh	:= yes
need_radphys	:= yes
need_linalg	:= yes
need_P1Diffusion := yes

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

# Put non template code here

libsrcs = Diffusion_DB.cc \
          P13T_Mesh_XYZ_pt.cc

# Get the library built, minding our various other targets too.

pretarg = deps

resolve_libs = \
	-L../mesh -lmesh_$(env) \
	-L../geom -lgeom_$(env) \
	-L../P1Diffusion -lP1Diffusion_$(env) \
	-L../radphys -lradphys_$(env) \
	-L../linalg -llinalg_$(env) \
	-L../util -lutil_$(env) \
	-L../spdf -lspdf_$(env) \
	-L../c4 -lc4_$(env) \
	-L../nml -lnml_$(env) \
	-L../ds++ -lds++_$(env)

include ../config/make.lib

deps:
ifneq ($(dodeps),no)
	(cd ../ds++;   $(MAKE) $(subnj))
	(cd ../nml;    $(MAKE) $(subnj) dodeps=no)
	(cd ../c4;     $(MAKE) $(subnj) dodeps=no)
	(cd ../spdf;   $(MAKE) $(subnj) dodeps=no)
	(cd ../util;   $(MAKE) $(subnj) dodeps=no)
	(cd ../linalg; $(MAKE) $(subnj) dodeps=no)
	(cd ../radphys; $(MAKE) $(subnj) dodeps=no)
	(cd ../P1Diffusion; $(MAKE) $(subnj) dodeps=no)
	(cd ../geom;   $(MAKE) $(subnj) dodeps=no)
	(cd ../mesh;   $(MAKE) $(subnj) dodeps=no)
endif

# Rules for building docs

docsrcs	= $(wildcard *.hh *.cc)

include ../config/make.doc

# Go pick up the standard signature for makefiles.

testclean:
	(cd test;     $(MAKE) clean)
	(cd testP13T; $(MAKE) clean)

aux_clean := testclean

cc_files := $(libsrcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/3T/makefile
###############################################################################
