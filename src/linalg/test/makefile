###############################################################################
# draco/src/linalg/test/makefile
# Dave Nystrom
# 14 January 1997
###############################################################################
# @> Rules for constructing libtstlinalg.a.
###############################################################################

pkg := tstlinalg

need_linalg   = yes
need_util     = yes
need_namelist = yes
need_spdf     = yes
need_c4       = yes
need_pcg      = yes
need_lapack   = yes
need_libm     = yes

SRC_TOP = ../..

###############################################################################
# Go get the standard includeable makefile for DRACO components.

include $(SRC_TOP)/config/makefile

export CXX CXXFLAGS PRST_CXXPP

###############################################################################
# Put non template code here

libsrcs = \
	tstpcg_DB.cc

###############################################################################
# Template files.  When is this really used.

tmplsrcs = 

###############################################################################
# Get the library built, minding our various other targets too.

pretarg = deps
auxtarg = $(TEST_PROGS)
resolve_libs = 

include $(SRC_TOP)/config/make.lib

###############################################################################
# Other targets.

deps:
ifneq ($(dodeps),no)
	(cd $(SRC_TOP)/ds++;	$(MAKE))
	(cd $(SRC_TOP)/c4;	$(MAKE) "dodeps=no" "do=")
	(cd $(SRC_TOP)/spdf;	$(MAKE) "dodeps=no" "do=")
	(cd $(SRC_TOP)/nml;	$(MAKE) "dodeps=no")
	(cd $(SRC_TOP)/util;	$(MAKE) "dodeps=no")
	(cd $(SRC_TOP)/linalg;	$(MAKE) "dodeps=no")
endif

###############################################################################
# Perstalyser targets.

tstpcg_DB_prs.cc:	tstpcg_DB.hh tstpcg_DB.prs .nml_tstpcg.hh

###############################################################################
# Now, build test programs.

tstpretarg = deps

tstsrcs = \
	TstPCG.cc

include $(SRC_TOP)/config/make.tests

###############################################################################
# Rules for building docs

docsrcs	= $(wildcard *.h *.cc)

include $(SRC_TOP)/config/make.doc

###############################################################################
# Go pick up the standard signature for makefiles.

cc_files := $(tstsrcs) $(libsrcs)

include $(SRC_TOP)/config/make.end

###############################################################################
# 			End of draco/src/linalg/test/makefile
###############################################################################
