###############################################################################
# draco/src/pydraco/makefile
# Geoffrey Furnish
# 8 January 1998
###############################################################################
# @> Rules for building the Draco/Python exectuable program.
###############################################################################

need_python	:= yes

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

srcs	= pydraco.cc

objs	= $(srcs:%.cc=$(env)/%.o)

default: pydraco_$(env)

pydraco_$(env) : deps $(objs) \
	$(dep_libds++)
	$(C++) $(C++FLAGS) -o $@ $(objs) $(LINK_LIBS)

deps:
ifneq ($(dodeps),no)
	(cd ../ds++; $(MAKE))
endif

# Go pick up the standard signature for makefiles.

cc_files := $(srcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/pydraco/makefile
###############################################################################
