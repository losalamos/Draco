###############################################################################
# draco/src/radphys/test/makefile
# Geoffrey Furnish
# 29 September 1997
###############################################################################
# @> Rules for building the baby solver test executable.
###############################################################################

SRC_TOP = ../..

need_units	:= yes
need_radphys	:= yes

# Go get the standard includeable makefile for Draco components.

include ../../config/makefile

srcs	= \
	testRadPhys.cc \
	testRadPhys_pt.cc \
	vectorops.cc

ifdef preinstantiate
# Add in files which perform explicit template instantiation.

pt_files = $(wildcard *.pt)
ptcc_files = $(pt_files:%.pt=.pt_%_$(env).cc)

endif

objs	= $(srcs:%.cc=$(env)/%.o) \
		$(ptcc_files:.%_$(env).cc=$(env)/%.o)

default: main_$(env)

main_$(env) : deps $(objs) \
	$(dep_libradphys) $(dep_libunits) $(dep_libds++)
	$(C++) $(C++FLAGS) -o $@ $(objs) $(LINK_LIBS) -lm

deps:
ifneq ($(dodeps),no)
	(cd ../../units;     $(MAKE) $(subnj))
	(cd ../../radphys;   $(MAKE) $(subnj) dodeps=no)
endif

# $(dep_libradphys):
$(dep_libds++):

# Go pick up the standard signature for makefiles.

aux_clean := exe_clean

exe_clean:
	$(RM) main_$(env)

cc_files := $(srcs)

include $(SRC_TOP)/config/make.end

###############################################################################
# 			End of draco/src/radphys/test/makefile
###############################################################################
