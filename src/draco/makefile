###############################################################################
# draco/src/draco/makefile
# Geoffrey Furnish
# 26 June 1997
###############################################################################
# @> Rules for building the Draco exectuable program.
###############################################################################

need_c4		:= yes
need_pooma	:= yes
need_mesh	:= yes

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

srcs	= main.cc draco.cc draconml.cc

ifdef preinstantiate
# Add in files which perform explicit template instantiation.

pt_files = $(wildcard *.pt)
ptcc_files = $(pt_files:%.pt=.pt_%_$(env).cc)

endif

objs	= $(srcs:%.cc=$(env)/%.o) \
		$(ptcc_files:.%_$(env).cc=$(env)/%.o)

#default: deps draco_$(env)
default: draco_$(env)

draco_$(env) : deps $(objs) \
	$(dep_libmesh) $(dep_libgeom) \
	$(dep_libnml) $(dep_libc4) $(dep_libds++)
	$(C++) $(C++FLAGS) -o $@ $(objs) $(LINK_LIBS)

deps:
ifneq ($(dodeps),no)
	(cd ../ds++; $(MAKE))
	(cd ../nml;  $(MAKE) dodeps=no)
	(cd ../c4;   $(MAKE) dodeps=no)
	(cd ../geom; $(MAKE) dodeps=no)
	(cd ../mesh; $(MAKE) dodeps=no)
endif

# Go pick up the standard signature for makefiles.

cc_files := $(srcs)

include ../config/make.end

###############################################################################
# 			End of draco/src/draco/makefile
###############################################################################
