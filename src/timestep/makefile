###############################################################################
# draco/src/timestep/makefile
# Geoffrey Furnish
# 4 December 1997
###############################################################################
# @> Rules for constructing the timestep support library for Draco (libtimestep.a).
###############################################################################

pkg := timestep

need_mesh	:= yes
need_linalg	:= yes

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

# Put non template code here

libsrcs = ts_manager.cc \
	  fixed_ts_advisor.cc \
	  ratio_ts_advisor.cc \
	  target_ts_advisor.cc \
          ts_advisor.cc

# Get the library built, minding our various other targets too.

pretarg = deps

resolve_libs = -L../ds++ -lds++_$(env)

include ../config/make.lib

deps:
ifneq ($(dodeps),no)
	(cd ../ds++;   $(MAKE) $(subnj))
endif

# Rules for building docs

docsrcs	= $(wildcard *.hh *.cc)

include ../config/make.doc

# Go pick up the standard signature for makefiles.

cc_files := $(libsrcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/timestep/makefile
###############################################################################
