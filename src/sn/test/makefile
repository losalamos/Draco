###############################################################################
# draco/src/sn/test/makefile
# Scott Turner - based on Geoff Furnish's makefiles
# 17 April 1998
###############################################################################
# @> Rules for building the solver test executable.
###############################################################################

SRC_TOP = ../..

need_sn := yes

# Go get the standard includeable makefile for Draco components.

include ../../config/makefile

srcs	= snpp.cc

objs	= $(srcs:%.cc=$(env)/%.o)

default: sn_$(env)

sn_$(env) : deps $(objs)
	$(C++) $(C++FLAGS) -o $@ $(objs) $(LINK_LIBS) -lm

deps:
ifneq ($(dodeps),no)
#	(cd ../../ds++;   $(MAKE) $(subnj))
	(cd ../../sn;     $(MAKE) $(subnj))
endif

# Go pick up the standard signature for makefiles.

aux_clean := exe_clean

exe_clean:
	$(RM) sn_$(env)

cc_files := $(srcs)

include $(SRC_TOP)/config/make.end

###############################################################################
# 			End of draco/src/sn/makefile
###############################################################################
