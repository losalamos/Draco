###############################################################################
# draco/src/3T/test/makefile
# Geoffrey Furnish
# 29 September 1997
###############################################################################
# @> Rules for building the baby solver test executable.
###############################################################################

SRC_TOP = ../..

need_linalg     := yes
need_c4		:= yes
need_namelist	:= yes
need_spdf	:= yes
need_util	:= yes
#need_pooma	:= yes
need_mesh	:= yes
need_3T		:= yes

# Go get the standard includeable makefile for Draco components.

include ../../config/makefile

srcs	= \
	main.cc \
	Quad_Params.cc \
	Trig_Params.cc \
	Test_Prob.cc \
	Run_DB.cc

ifdef preinstantiate
# Add in files which perform explicit template instantiation.

pt_files = $(wildcard *.pt)
ptcc_files = $(pt_files:%.pt=.pt_%_$(env).cc)

endif

objs	= $(srcs:%.cc=$(env)/%.o) \
		$(ptcc_files:.%_$(env).cc=$(env)/%.o)

default: main_$(env)

main_$(env) : deps $(objs) \
	$(dep_lib3T) $(dep_libmesh) $(dep_libgeom) $(dep_liblinalg) \
	$(dep_libutil) $(dep_libspdf) $(dep_libnml) $(dep_libc4) $(dep_libds++)
	$(C++) $(C++FLAGS) -o $@ $(objs) $(LINK_LIBS) -lm

deps:
ifneq ($(dodeps),no)
	(cd ../../3T;   $(MAKE) $(subnj))
	(cd ../../spdf;   $(MAKE) $(subnj) dodeps=no)
	(cd ../../util;   $(MAKE) $(subnj) dodeps=no)
endif

$(dep_lib3T):
$(dep_libmesh):
$(dep_libgeom):
$(dep_liblinalg):
$(dep_libutil):
$(dep_libspdf):
$(dep_libnml):
$(dep_libc4):
$(dep_libds++):

# Go pick up the standard signature for makefiles.

aux_clean := exe_clean

exe_clean:
	$(RM) main_$(env)

cc_files := $(srcs)

include $(SRC_TOP)/config/make.end

###############################################################################
# 			End of draco/src/3T/test/makefile
###############################################################################
