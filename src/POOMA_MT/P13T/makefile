###############################################################################
# draco/src/POOMA_MT/P13T/makefile
# Randy M. Roberts
# 29 September 1997
###############################################################################
# @> Rules for building the P13T instantiations
###############################################################################

pkg := PoomaP13T

SRC_TOP = ../..

need_pooma := yes
need_c4	   := yes
need_3T    := yes
need_timestep := yes
need_units	:= yes
need_radphys := yes
need_linalg	:= yes
need_P1Diffusion := yes

strict := no

# Go get the standard includeable makefile for Draco components.

include $(SRC_TOP)/config/makefile

include $(SRC_TOP)/config/.pooma.$(env)

libsrcs	= \
	PreCondP1Diff_pt.cc \
	MatVecP1Diff_pt.cc \
	SolverP1Diff_pt.cc \
	radphys_pt.cc \
	P1Diffusion_pt.cc \
	MatrixP1Diff_pt.cc \
	P13T_PoomaMesh_XYZ_pt.cc

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

# Get the library built, minding our various other targets too.

pretarg = deps

resolve_libs = \
	-L$(SRC_TOP)/mesh -lmesh_$(env) \
	-L$(SRC_TOP)/geom -lgeom_$(env) \
	-L$(SRC_TOP)/P1Diffusion -lP1Diffusion_$(env) \
	-L$(SRC_TOP)/radphys -lradphys_$(env) \
	-L$(SRC_TOP)/linalg -llinalg_$(env) \
	-L$(SRC_TOP)/util -lutil_$(env) \
	-L$(SRC_TOP)/spdf -lspdf_$(env) \
	-L$(SRC_TOP)/c4 -lc4_$(env) \
	-L$(SRC_TOP)/nml -lnml_$(env) \
	-L$(SRC_TOP)/ds++ -lds++_$(env)

include $(SRC_TOP)/config/make.lib

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
# 			End of draco/src/3T/P13T/makefile
###############################################################################
