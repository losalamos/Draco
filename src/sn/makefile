###############################################################################
# draco/src/sn/makefile
# Scott Turner - based on Geoff Furnish's makefiles
# 19 February 1998
###############################################################################
# @> Rules for constructing the sn support library for Draco (libsn.a).
###############################################################################

pkg := sn

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

# Put non template code here

libsrcs = inner.cc sweep3d.cc

# Get the library built, minding our various other targets too.

include ../config/make.lib

# Rules for building docs

docsrcs	= $(wildcard *.hh *.cc)

include ../config/make.doc

# Go pick up the standard signature for makefiles.

cc_files := $(libsrcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/sn/makefile
###############################################################################
