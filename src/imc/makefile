###############################################################################
# draco/src/imc/makefile 
# Thomas M. Evans
# Tue Mar 10 11:44:16 1998
###############################################################################
# @> Rules for constructing the IMC support library
###############################################################################

pkg := imc

need_c4  := yes
need_rng := yes

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

# Code sources

# tempsrc = template code
# nonsrc  = non template code
# libsrcs = code to be compiled, if using instantialyzer only include
#           nonsrc

nonsrc  = Layout.cc OS_Builder.cc OS_Mesh.cc OS_Interface.cc Coord_sys.cc

tempsrc = 

libsrcs = $(nonsrc) $(tempsrc)

# Get the library built, minding our various other targets too.

pretarg = deps
resolve_libs = $(dep_libds++) $(dep_libc4) $(dep_librng)

include ../config/make.lib

deps:
ifneq ($(dodeps),no)
	(cd ../ds++; $(MAKE))
	(cd ../c4;   $(MAKE) dodeps=no)
	(cd ../rng;  $(MAKE) dodeps=no)
endif

# Rules for building docs

docsrcs = $(wildcard *.hh *.cc)

include ../config/make.doc

# Go pick up the standard signature for makefiles.

cc_files := $(libsrcs)

include ../config/make.end

###############################################################################
#			 end of draco/src/imc/makefile
###############################################################################
