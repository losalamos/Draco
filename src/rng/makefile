###############################################################################
# draco/src/rng/makefile 
# Thomas M. Evans
# Wed Apr 29 11:02:46 1998
###############################################################################
# @> Rules for constructing rng support library
###############################################################################

pkg := rng

need_sprng := yes

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

# Code sources

# tempsrc = template code
# nonsrc  = non template code
# libsrcs = code to be compiled, if using instantialyzer only include
#           nonsrc

nonsrc  = SMrng.cc Rnd_Control.cc

tempsrc =

libsrcs = $(nonsrc) $(tempsrc)

# Get the library built, minding our various other targets too.

pretarg = deps
resolve_libs = $(dep_libds++)

include ../config/make.lib

deps:
ifneq ($(dodeps),no)
	(cd ../ds++; $(MAKE))
endif

# Rules for building docs

docsrcs = $(wildcard *.hh *.cc)

include ../config/make.doc

# Go pick up the standard signature for makefiles.

cc_files := $(libsrcs)

include ../config/make.end

###############################################################################
#			 end of draco/src/rng/makefile
###############################################################################
