###############################################################################
# draco/src/timestep/test/makefile
# John McGhee
# 3 Jun 1998
###############################################################################
# @> Rules for building the timestep test executable.
###############################################################################

pkg	:= timestep

SRC_TOP = ../..

#need_linalg    := yes
#need_c4	:= yes
#need_namelist	:= yes
#need_spdf	:= yes
#need_util	:= yes
#need_pooma	:= yes
#need_mesh	:= yes
need_timestep	:= yes

# Go get the standard includeable makefile for Draco components.

include ../../config/makefile

srcs	= \
	main.cc \
	test_timestep.cc \
	test_timestep_pt.cc \
	dummy_package.cc   \
	test_utils.cc

ifdef preinstantiate
# Add in files which perform explicit template instantiation.

pt_files = $(wildcard *.pt)
ptcc_files = $(pt_files:%.pt=.pt_%_$(env).cc)

endif

objs	= $(srcs:%.cc=$(env)/%.o) \
		$(ptcc_files:.%_$(env).cc=$(env)/%.o)

default: main_$(env)

main_$(env) : deps $(objs) \
	$(dep_libtimestep) $(dep_libmesh) $(dep_libgeom) $(dep_liblinalg) \
	$(dep_libutil) $(dep_libspdf) $(dep_libnml) $(dep_libc4) $(dep_libds++)
	$(C++) $(C++FLAGS) -o $@ $(objs) $(LINK_LIBS) -lm

deps:
ifneq ($(dodeps),no)
	(cd ../../timestep;   $(MAKE) $(subnj))
	(cd ../../spdf;   $(MAKE) $(subnj) dodeps=no)
	(cd ../../util;   $(MAKE) $(subnj) dodeps=no)
endif

$(dep_libtimestep):
$(dep_libmesh):
$(dep_libgeom):
$(dep_liblinalg):
$(dep_libutil):
$(dep_libspdf):
$(dep_libnml):
$(dep_libc4):
$(dep_libds++):

# Get DejaGnu support

include ../../config/make.deja

# Go pick up the standard signature for makefiles.

aux_clean := exe_clean

exe_clean:
	$(RM) main_$(env)

cc_files := $(srcs)

include $(SRC_TOP)/config/make.end

###############################################################################
# 			End of draco/src/timestep/test/makefile
###############################################################################
