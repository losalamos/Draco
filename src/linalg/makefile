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
	-L../util -lutil_$(env) \
	-L../spdf -lspdf_$(env) \
	-L../c4 -lc4_$(env) \
	-L../nml -lnml_$(env) \
	-L../ds++ -lds++_$(env)

include ../config/make.lib

deps:
ifneq ($(dodeps),no)
	(cd ../ds++; $(MAKE))
	(cd ../nml;  $(MAKE) $(subnj) dodeps=no)
	(cd ../c4;   $(MAKE) $(subnj) dodeps=no)
	(cd ../spdf; $(MAKE) $(subnj) dodeps=no)
	(cd ../util; $(MAKE) $(subnj) dodeps=no)
endif

###############################################################################
# The test programs.

tstsrcs	:= $(wildcard t*.cc)
tstpretarg := $(pretarg)

include ../config/make.tests

###############################################################################
# Testing out some dejagnu support

deja:
	(cd testsuites; $(MAKE) links)
	runtest --tool=linalg --srcdir=./testsuites ENV=$(env) DBC=$(DBC)

###############################################################################
# Perstalyser targets.

pcg_DB_prs.cc: pcg_DB.hh pcg_DB.prs .nml_pcg.hh

###############################################################################
# Rules for building docs

docsrcs	= $(wildcard *.h *.cc)

include ../config/make.doc

###############################################################################
# Go pick up the standard signature for makefiles.

testclean:
	(cd test; $(MAKE) clean)

aux_clean := testclean

cc_files := $(libsrcs) $(tstsrcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/linalg/makefile
###############################################################################
