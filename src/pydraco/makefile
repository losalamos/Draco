###############################################################################
# draco/src/pydraco/makefile
# Geoffrey Furnish
# 8 January 1998
###############################################################################
# @> Rules for building the Draco/Python exectuable program.
###############################################################################

need_python	:= yes
need_util	:= yes

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

srcs	= pydraco.cc adfilemodule.cc #cxxmodule.cc

objs	= $(srcs:%.cc=$(env)/%.o)

default: deps pydraco_$(env) adfilemodule.so

pydraco_$(env) : $(objs) \
	$(dep_libds++) $(dep_libc4) $(dep_libspdf) $(dep_libutil)
	$(C++) $(C++FLAGS) -o $@ $(objs) $(LINK_LIBS)

adfilemodule.so: $(env)/adfilemodule.o
	$(LDSHARED) -o $@ $(env)/adfilemodule.o

deps:
ifneq ($(dodeps),no)
	(cd ../ds++; $(MAKE))
	(cd ../c4;   $(MAKE))
	(cd ../spdf; $(MAKE))
	(cd ../util; $(MAKE))
endif

# Go pick up the standard signature for makefiles.

cc_files := $(srcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/pydraco/makefile
###############################################################################
