###############################################################################
# draco/src/linalg/makefile
# Dave Nystrom
# 13 August 1996
###############################################################################
# @> Rules for constructing liblinalg.a.
###############################################################################

pkg := linalg

need_util = yes

###############################################################################
# Go get the standard includeable makefile for DRACO components.

include ../config/makefile

export CXX CXXFLAGS PRST_CXXPP

###############################################################################
# Put non template code here

libsrcs = \
	PCG_Subroutines.cc \
	pcg_DB.cc

###############################################################################
# Template files.

tmplsrcs = 

###############################################################################
# Get the library built, minding our various other targets too.

pretarg = deps
auxtarg = 
resolve_libs = \
	../util/libutil_$(env).$(lib_suffix) \
	../spdf/libspdf_$(env).$(lib_suffix) \
	../c4/libc4_$(env).$(lib_suffix) \
	../nml/libnml_$(env).$(lib_suffix) \
	../ds++/libds++_$(env).$(lib_suffix)

include ../config/make.lib

deps:
ifneq ($(dodeps),no)
	(cd ../ds++; $(MAKE))
	(cd ../nml;  $(MAKE) dodeps=no)
	(cd ../c4;   $(MAKE) dodeps=no)
	(cd ../spdf; $(MAKE) dodeps=no)
	(cd ../util; $(MAKE) dodeps=no)
endif

###############################################################################
# The test programs.

tstsrcs	:= $(wildcard t*.cc)
tstpretarg := $(pretarg)

include ../config/make.tests

###############################################################################
# Perstalyser targets.

pcg_DB_prs.cc: pcg_DB.hh pcg_DB.prs .nml_pcg.hh

###############################################################################
# Rules for building docs

docsrcs	= $(wildcard *.h *.cc)

include ../config/make.doc

###############################################################################
# Go pick up the standard signature for makefiles.

cc_files := $(libsrcs) $(tstsrcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/linalg/makefile
###############################################################################
