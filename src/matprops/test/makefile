###############################################################################
# draco/src/matprops/test/makefile
# Geoffrey Furnish
# 29 September 1997
###############################################################################
# @> Rules for building the baby solver test executable.
###############################################################################

SRC_TOP = ../..

need_linalg     := no
need_c4		:= no
need_namelist	:= no
need_spdf	:= no
need_util	:= no
need_pooma	:= no
need_mesh	:= no
need_units	:= yes
need_matprops	:= yes


# Go get the standard includeable makefile for Draco components.

include ../../config/makefile

# turn on full debugging
CXXFLAGS += -DDBC=7

srcs	= \
	main.cc \
	testBilin.cc \
	testMatProp.cc \
	testMatProp_pt.cc

ifdef preinstantiate
# Add in files which perform explicit template instantiation.

pt_files = $(wildcard *.pt)
ptcc_files = $(pt_files:%.pt=.pt_%_$(env).cc)

endif

objs	= $(srcs:%.cc=$(env)/%.o) \
		$(ptcc_files:.%_$(env).cc=$(env)/%.o)

default: main_$(env)

main_$(env) : deps $(objs) \
	$(dep_libmatprops) $(dep_libds++)
	$(C++) $(C++FLAGS) -o $@ $(objs) $(LINK_LIBS) -lm

deps:
ifneq ($(dodeps),no)
	(cd ../../matprops;   $(MAKE) $(subnj))
endif

$(dep_libmatprops):
$(dep_libds++):

# Go pick up the standard signature for makefiles.

aux_clean := exe_clean

exe_clean:
	$(RM) main_$(env)

cc_files := $(srcs)

include $(SRC_TOP)/config/make.end

###############################################################################
# 			End of draco/src/matprops/test/makefile
###############################################################################
