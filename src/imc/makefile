###############################################################################
# draco/src/imctest/makefile 
# Thomas M. Evans
# Tue Mar 10 11:44:16 1998
###############################################################################
# @> Rules for constructing the IMCTEST support library
###############################################################################

pkg := imctest

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

# Put non template code here

libsrcs = Layout.cc OS_Builder.cc OS_Mesh.cc OS_Parser.cc Random.cc \
          XYCoord_sys.cc XYZCoord_sys.cc

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
#			 end of draco/src/imctest/makefile
###############################################################################
