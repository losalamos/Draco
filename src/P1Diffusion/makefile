###############################################################################
# draco/src/P1Diffusion/makefile
# Randy M. Roberts
# 01 October 1998
###############################################################################
# @> Rules for constructing the P1Diffusion support library for Draco (libP1Diffusion.a).
###############################################################################

pkg := P1Diffusion

need_ds++	:= yes
need_linalg	:= yes
need_mesh	:= yes
need_diffusion  := yes

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

# Put non template code here

libsrcs = MatrixP1Diff_pt.cc \
	P1Diffusion_pt.cc \
	PreCondP1Diff_pt.cc \
	SolverP1Diff_pt.cc \
	MatVecP1Diff_pt.cc

# Get the library built, minding our various other targets too.

pretarg = deps

resolve_libs = -L../ds++ -lds++_$(env) \
               -L../linalg -llinalg_$(env) \
               -L../mesh -lmesh_$(env) \
               -L../diffusion -ldiffusion_$(env)

include ../config/make.lib

deps:
ifneq ($(dodeps),no)
	(cd ../ds++;    $(MAKE) $(subnj))
	(cd ../linalg;    $(MAKE) $(subnj))
	(cd ../mesh;    $(MAKE) $(subnj))
	(cd ../diffusion;    $(MAKE) $(subnj))
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
# 			End of draco/src/P1Diffusion/makefile
###############################################################################
