###############################################################################
# draco/src/3T/test/makefile
# Geoffrey Furnish
# 29 September 1997
###############################################################################
# @> Rules for building the baby solver test executable.
###############################################################################

SRC_TOP = ../..

need_pooma := yes
strict := no

# Go get the standard includeable makefile for Draco components.

include ../../config/makefile

include $(SRC_TOP)/config/.pooma.$(env)

srcs	= \
	Pooma_txyz.cc

# CXXFLAGS += --restrict --output_restrict --backend -LANG:restrict=on \
#	$(shell ./stripPoomaDefines.sh ${POOMA_ROOT} ${POOMA_ARCH}) \
#	-I${POOMA_ROOT}/src

CXXFLAGS += $(POOMA_INCLUDES) $(POOMA_DEFINES)

# LINK_LIBS += -L${POOMA_ROOT}/lib/${POOMA_ARCH} \
#	-lpooma_templates_basic -lpooma_templates_3D -lpooma -lfft

LINK_LIBS	+= $(POOMA_LFLAGS) $(POOMA_LIBS)

ifdef preinstantiate
# Add in files which perform explicit template instantiation.

pt_files = $(wildcard *.pt)
ptcc_files = $(pt_files:%.pt=.pt_%_$(env).cc)

endif

objs	= $(srcs:%.cc=$(env)/%.o) \
		$(ptcc_files:.%_$(env).cc=$(env)/%.o)

default: main

main: deps main_$(env)
	@echo 
	@echo "main target has been built"
	@echo 

main_$(env) : $(objs) \
	$(dep_libds++) \
	$(dep_libpooma)
	$(C++) $(C++FLAGS) \
		-o $@ $(objs) $(LINK_LIBS) -lm

deps:
ifneq ($(dodeps),no)

endif

$(dep_libds++):
$(dep_libpooma):

# Go pick up the standard signature for makefiles.

aux_clean := exe_clean

exe_clean:
	$(RM) main_$(env)
	$(RM) $(SRC_TOP)/config/.pooma.$(env)

cc_files := $(srcs)

include $(SRC_TOP)/config/make.end

###############################################################################
# 			End of draco/src/3T/test/makefile
###############################################################################
