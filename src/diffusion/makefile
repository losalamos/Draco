###############################################################################
# draco/src/diffusion/makefile
# Randy M. Roberts
# 01 October 1998
###############################################################################
# @> Rules for constructing the diffusion support library for Draco (libdiffusion.a).
###############################################################################

pkg := diffusion

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

# Put non template code here

libsrcs = Diffusion_DB.cc

# Get the library built, minding our various other targets too.

pretarg = deps

resolve_libs = -L../ds++ -lds++_$(env) \
               -L../linalg -llinalg_$(env) \
	       -L../nml -lnml_$(env) \
               -L../mesh -lmesh_$(env)

include ../config/make.lib

deps:
ifneq ($(dodeps),no)
	(cd ../nml;    $(MAKE) $(subnj))
endif

# Rules for building docs

docsrcs	= $(wildcard *.hh *.cc)

include ../config/make.doc

# Go pick up the standard signature for makefiles.

testclean:

aux_clean := testclean

cc_files := $(libsrcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/diffusion/makefile
###############################################################################
