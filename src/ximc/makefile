###############################################################################
# draco/src/imc/test/makefile 
# Thomas M. Evans
# Tue Mar 10 15:00:20 1998
###############################################################################
# @> Rules for building the testbuilder
###############################################################################

SRC_TOP = ..

need_imc        := yes
need_c4         := yes

# Go get the standard includeable makefile for Draco components.

include ../config/makefile

srcs = imctest.cc

ifdef preinstantiate
# Add in files which perform explicit template instantiation.

pt_files = $(wildcard *.pt)
ptcc_files = $(pt_files:%.pt=.pt_%_$(env).cc)

endif

objs    = $(srcs:%.cc=$(env)/%.o) \
                $(ptcc_files:.%_$(env).cc=$(env)/%.o)

default: main_$(env)

main_$(env) : deps $(objs) \
	$(dep_libimc) $(dep_libds++) $(dep_libc4)
	$(C++) $(C++FLAGS) -o $@ $(objs) $(LINK_LIBS) -lm 

deps:
ifneq ($(dodeps),no)
	(cd ../imc;   $(MAKE) $(subnj))
endif

$(dep_libimc):
$(dep_libds++):
$(dep_libc4):

# Go pick up the standard signature for makefiles.

aux_clean := exe_clean

exe_clean:
	$(RM) main_$(env)

cc_files := $(srcs)

include $(SRC_TOP)/config/make.end


###############################################################################
#			 end of draco/src/imc/test/makefile
###############################################################################
