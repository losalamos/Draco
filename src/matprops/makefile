###############################################################################
# draco/src/matprops/makefile
# Randy M. Roberts
# 7 April 1998
###############################################################################
# @> Rules for constructing the matprops support library for Draco (libmatprops.a).
###############################################################################

pkg := matprops

need_radphys	:= yes

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

# turn on full debugging
CXXFLAGS += -DDBC=7

# Put non template code here

libsrcs = BilinearInterpTable.cc BilinearInterpGrid.cc \
          FifiParser.cc FifiMatPropsReader.cc InterpedMaterialProps.cc \
          TempMapper_pt.cc InterpedMaterialProps_pt.cc  \
	  MarshakMaterialProps_pt.cc

# Get the library built, minding our various other targets too.

pretarg = deps

resolve_libs = \
	-L../radphys -lradphys_$(env) \
	-L../units -lunits_$(env) \
	-L../ds++ -lds++_$(env)

include ../config/make.lib

deps:
ifneq ($(dodeps),no)
	(cd ../ds++;    $(MAKE) $(subnj))
	(cd ../units;   $(MAKE) $(subnj) dodeps=no)
	(cd ../radphys; $(MAKE) $(subnj) dodeps=no)
endif

$(dep_libradphys):

# Rules for building docs

docsrcs	= $(wildcard *.hh *.cc)

include ../config/make.doc

# Go pick up the standard signature for makefiles.

cc_files := $(libsrcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/matprops/makefile
###############################################################################
