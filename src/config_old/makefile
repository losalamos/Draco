###############################################################################
# draco/src/config/makefile
# Geoffrey Furnish
# 25 July 1997
###############################################################################
# @> This "makefile" is included by the makefiles in each component
# @> directory.  The intention is that this makefile will hold all the
# @> definitions required to configure compilers and tools for building
# @> the Draco components.  The makefile for each component will only need
# @> to include this one and then do its own component specific stuff
# @> without having to worry about how to set up all the compile and link
# @> options. 
###############################################################################

# Get support for parallel makes configured.

ifndef nj
nj	:= 1
endif

subnj	:= -j $(nj)

# Now setup the suffixes we are prepared to work with.

.SUFFIXES:
.SUFFIXES: .cpp .cc .o

# Find out where we are running.

uname = $(shell uname)
host  = $(shell hostname)

# Set up some defaults so that lossage will be obvious if we blow
# something. 

CC = bogus
C++ = bogus
C++FLAGS = bogus
C++HDR_DIR = bogus
PRIVATE_LIB_DIR	= $(HOME)/lib

# Fix HOSTTYPE macro on Linux if some crazy shell emulator breaks it.
# Lemacs 19.10, for instance, sets it to i386-linux.

ifneq (,$(findstring i386,${HOSTTYPE}))
HOSTTYPE := i386
endif

# Ask not to inherit certain variables.

unexport need_gtslib need_gf++ need_tcltk need_namelist need_plplot need_x11 need_python

#-----------------------------------------------------------------------------#
# Configure for Physical Dynamics C++ Namelist Utility.

NML_LIB_DIR	= $(SRC_TOP)/nml

ifeq ($(need_tcltk),yes)
NML_LIBS	= -L$(NML_LIB_DIR) -lnmltk_$(env)
else
NML_LIBS	= -L$(NML_LIB_DIR) -lnml_$(env)
endif

#-----------------------------------------------------------------------------#
# Configure for C4 (Cannonical Classes for Concurrency Control)

C4_LIB		= -L$(SRC_TOP)/c4 -lc4_$(env)

#-----------------------------------------------------------------------------#
# Configure for Draco development.

ifndef SRC_TOP
SRC_TOP	= ..
endif

DS++_LIB	= -L$(SRC_TOP)/ds++     -lds++_$(env)

SPDF_LIB	= -L$(SRC_TOP)/spdf     -lspdf_$(env)
UTIL_LIB	= -L$(SRC_TOP)/util     -lutil_$(env)
LINALG_LIB	= -L$(SRC_TOP)/linalg   -llinalg_$(env)
GEOM_LIB	= -L$(SRC_TOP)/geom     -lgeom_$(env)
MESH_LIB	= -L$(SRC_TOP)/mesh     -lmesh_$(env)
3T_LIB		= -L$(SRC_TOP)/3T       -l3T_$(env)
TIMESTEP_LIB 	= -L$(SRC_TOP)/timestep -ltimestep_$(env)
UNITS_LIB	= -L$(SRC_TOP)/units    -lunits_$(env)
RADPHYS_LIB	= -L$(SRC_TOP)/radphys  -lradphys_$(env)
MATPROPS_LIB	= -L$(SRC_TOP)/matprops -lmatprops_$(env)
SN_LIB		= -L$(SRC_TOP)/sn       -lsn_$(env)
IMC_LIB         = -L$(SRC_TOP)/imc      -limc_$(env)
PYCXX_LIB	= -L$(SRC_TOP)/pycxx    -lpycxx_$(env)
RNG_LIB         = -L$(SRC_TOP)/rng      -lrng_$(env)
P1DIFFUSION_LIB  = -L$(SRC_TOP)/P1Diffusion -lP1Diffusion_$(env)

ifeq ($(need_radphys),yes)
need_units	:= yes
endif

ifeq ($(need_util),yes)
need_spdf	:= yes
need_c4		:= yes
endif

ifeq ($(need_mesh),yes)
need_geom	:= yes
need_c4		:= yes
need_namelist	:= yes
endif

ifeq ($(need_linalg),yes)
need_pcg	:= yes
endif

ifeq ($(need_pcg),yes)
need_lapack	:= yes
endif

ifeq ($(need_imc),yes)
need_rng        := yes
endif

ifeq ($(need_rng),yes)
need_sprng      := yes
endif

#-----------------------------------------------------------------------------#
# Set up some tools.

MAKEDEPEND	:= gcc -MM -U__GNUC__

NMLGEN		:= $(SRC_TOP)/../tools/nmlgen

PERSTALYZER	:= $(SRC_TOP)/../tools/perstalyzer

INSTANTIALYZER	:= $(SRC_TOP)/../tools/instantialyzer

AUTODOC		:= $(SRC_TOP)/../tools/autodoc

FIXDEPS		:= $(SRC_TOP)/../tools/fixdeps

ifeq ($(prstdbg),on)
PERSTALYZER	+= -v
endif

# Turn on debugging by default, so core's can be debugged without the
# frustration of having to recompile to get debugging going.  During
# development, we take the philosophy the code should be easily
# debuggable by default, and you should have to take action if you
# /don't/ want debugging.

ifeq ($(debug),off)
DEBUG_FLAG	= 
else
DEBUG_FLAG	= -g
endif

# By default we leave assertions enabled, but we allow the user to
# turn them off he really wants to (by saying "make assert=no").

ifeq ($(assert),no)
ASSERT_FLAG	= -DNOASSERT
else
ASSERT_FLAG	=
endif

###############################################################################
# Now we perform any adjustments needed on particular platofrms.
###############################################################################

###############################################################################
# Go pick up system specific info.
# These files should define the following:
#	CC	C compiler
#	CXX	C++ compiler
#	env	Name of compilation environemnt (keep it short!)
#	MAKEDEPEND	Probably $(CXX) -M, or suchlike.
# and some other, less important, stuff which I will try to document
# one of these days.
#
# Note that although this first fan out is done on the basis of the OS
# name as reported by uname, it is also perfectly okay to fan out
# further within these files, perhpas on the basis of hostname, etc,
# since different physical machines may have different licensed
# compiler products, etc.

oscfgfile = $(SRC_TOP)/config/$(uname).cf

ifeq ($(wildcard $(oscfgfile)),$(oscfgfile))
include $(oscfgfile)
endif

hstcfgfile = $(SRC_TOP)/config/$(host).cf

ifeq ($(wildcard $(hstcfgfile)),$(hstcfgfile))
include $(hstcfgfile)
endif

###############################################################################
# Now setup the main compilation control variables.

# For now, just put all conceivable -I's and so forth on the compile
# line so that we can hopefully minimize the Cfront template
# instantiation trauma.

# This used to be in CXX_INC_FLAGS: -I$(GTS_SRC_TOP)/ds++ 

CXX_INC_FLAGS = -I$(SRC_TOP) $(POOMA_FLAGS) $(PY_INC_FLAG)

CXXFLAGS = $(SYS_FLAGS) $(CXX_INC_FLAGS) $(DEBUG_FLAG) $(XCPT_FLAG) \
	 $(ASSERT_FLAG) $(PYTHON_FLAG)

# Need these because I have some cases of explicit rules in my various
# makefiles, written back when I was using gmake before the
# conversion from C++ to CXX, etc.

C++ = $(CXX)
C++FLAGS = $(CXXFLAGS)

#-----------------------------------------------------------------------------#
# Now construct the macro used for linking programs.  Note that it is
# imperative that this be constructed in the proper order.  If x has a 
# dependency on y, then x must occur before y in the link line!

LINK_LIBS	:=
ifdef SYS_LINK_FLAGS
LINK_LIBS	+= $(SYS_LINK_FLAGS)
endif
ifeq ($(need_sn),yes)
LINK_LIBS	+= $(SN_LIB)
endif
ifeq ($(need_imc),yes)
LINK_LIBS       += $(IMC_LIB)
endif
ifeq ($(need_rng),yes)
LINK_LIBS       += $(RNG_LIB)
endif
ifeq ($(need_P1Diffusion),yes)
LINK_LIBS       += $(P1DIFFUSION_LIB)
endif
ifeq ($(need_matprops),yes)
LINK_LIBS	+= $(MATPROPS_LIB)
endif
ifeq ($(need_3T),yes)
LINK_LIBS	+= $(3T_LIB)
endif
ifeq ($(need_radphys),yes)
LINK_LIBS	+= $(RADPHYS_LIB)
endif
ifeq ($(need_units),yes)
LINK_LIBS	+= $(UNITS_LIB)
endif
ifeq ($(need_timestep),yes)
LINK_LIBS	+= $(TIMESTEP_LIB)
endif
ifeq ($(need_mesh),yes)
LINK_LIBS	+= $(MESH_LIB)
endif
ifeq ($(need_geom),yes)
LINK_LIBS	+= $(GEOM_LIB)
endif
ifeq ($(need_linalg),yes)
LINK_LIBS	+= $(LINALG_LIB)
endif
ifeq ($(need_pycxx),yes)
LINK_LIBS	+= $(PYCXX_LIB)
endif
ifeq ($(need_util),yes)
LINK_LIBS	+= $(UTIL_LIB)
endif
ifeq ($(need_spdf),yes)
LINK_LIBS	+= $(SPDF_LIB)
endif
ifeq ($(need_namelist),yes)
LINK_LIBS	+= $(NML_LIBS)
endif
ifeq ($(need_c4),yes)
LINK_LIBS	+= $(C4_LIB)
endif
LINK_LIBS	+= $(DS++_LIB)
ifeq ($(need_pcg),yes)
LINK_LIBS	+= $(PCG_LIB)
endif
ifeq ($(need_lapack),yes)
LINK_LIBS	+= $(LAPACK_LIB)
endif
ifeq ($(need_pooma),yes)
#LINK_LIBS	+= $(POOMA_LIB)
LINK_LIBS	+= $(POOMA_LFLAGS) $(POOMA_LIBS)
endif
ifeq ($(need_python),yes)
LINK_LIBS	+= $(PYLIBS)
need_libm	:= yes
endif
ifdef MPI_LIB
LINK_LIBS	+= $(MPI_LIB)
endif
ifdef SYS_LIBS
LINK_LIBS	+= $(SYS_LIBS)
endif
ifeq ($(need_libm),yes)
LINK_LIBS	+= -lm
endif

#-----------------------------------------------------------------------------#
# Now construct some macros for use in writing dependencies upon libraries.

ifdef make_shared_lib
lib_suffix	= $(shlib_suffix)
else
lib_suffix	= a
endif

dep_libds++	= $(SRC_TOP)/ds++/libds++_$(env).$(lib_suffix)
dep_libnml	= $(SRC_TOP)/nml/libnml_$(env).$(lib_suffix)
dep_libnmltk	= $(SRC_TOP)/nml/libnmltk_$(env).$(lib_suffix)
dep_libc4	= $(SRC_TOP)/c4/libc4_$(env).$(lib_suffix)
dep_libspdf	= $(SRC_TOP)/spdf/libspdf_$(env).$(lib_suffix)
dep_libutil	= $(SRC_TOP)/util/libutil_$(env).$(lib_suffix)
dep_libgeom	= $(SRC_TOP)/geom/libgeom_$(env).$(lib_suffix)
dep_libmesh	= $(SRC_TOP)/mesh/libmesh_$(env).$(lib_suffix)
dep_liblinalg	= $(SRC_TOP)/linalg/liblinalg_$(env).$(lib_suffix)
dep_lib3T	= $(SRC_TOP)/3T/lib3T_$(env).$(lib_suffix)
dep_libtimestep	= $(SRC_TOP)/timestep/libtimestep_$(env).$(lib_suffix)
dep_libunits	= $(SRC_TOP)/units/libunits_$(env).$(lib_suffix)
dep_libradphys	= $(SRC_TOP)/radphys/libradphys_$(env).$(lib_suffix)
dep_libmatprops	= $(SRC_TOP)/matprops/libmatprops_$(env).$(lib_suffix)
dep_libimc      = $(SRC_TOP)/imc/libimc_$(env).$(lib_suffix)
dep_librng      = $(SRC_TOP)/rng/librng_$(env).$(lib_suffix)
dep_libP1Diffusion = $(SRC_TOP)/P1Diffusion/libP1Diffusion_$(env).$(lib_suffix)

#dep_libpooma	= $(POOMA_LIB_DIR)/libpooma.a
dep_libpooma	= $(POOMA_ROOT)/lib/$(POOMA_ARCH)/libpooma.a

#-----------------------------------------------------------------------------#
# Now finally, some implicit rules which we will be using constantly,
# which depend critically on the value of $(env).

$(env)/%.o : %.cc
	$(CXX) -c $(CXXFLAGS) $(CXX_DEPEND) $*.cc -o $@
ifdef CXX2
	$(CXX2)
endif

.pt_%_$(env).cc : %.pt
	$(INSTANTIALYZER) -o $@ $*.pt

$(env)/%.o : .%_$(env).cc
	$(CXX) -c $(CXXFLAGS) $(CXX_DEPEND) .$*_$(env).cc -o $@
ifdef CXX2
	$(CXX2)
endif

.nml_%.hh .nml_%.cc : %.nml
	$(NMLGEN) $*.nml

%_prs.cc : %.prs
	$(PERSTALYZER) -c $*

# Use GNU style auto dependencies.  Read info on GNU make to
# understand.  Note that we want the dependency generation to reflect
# an honest status code, and the method shown in the GNU make manual
# does not do this.  The problem is that the exit status of a pipeline
# is the exit status of the last command in the pipeline.  In order to
# see if the dependency analysis phase (MAKEDEPEND) had problems, we
# have to save it's exit status, and return that.  In between, we
# conditionally generate the requested dependency file.

ifndef kcc_deps
ifndef sgi_deps
.d.%_$(env):	%.cc
	@$(SHELL) -c 'echo "Dependalyzing $@"; \
	$(MAKEDEPEND) $(CXX_DEPEND_FLAG) $< >.x.$@; \
	rc=$$?; \
	test $$rc = "0" && \
		(cat .x.$@ | sed '\''s/$*.o/$(env)\/& $@/g'\'' > $@); \
	rm .x.$@; \
	exit $$rc'

# A rule for dependalyzing template instantiation files.

.d.pt_%: .pt_%.cc
	@$(SHELL) -c 'echo "Dependalyzing $@"; \
	$(MAKEDEPEND) $(CXX_DEPEND_FLAG) $< > .x.$@; \
	rc=$$?; \
	test $$rc = "0" && \
		(perl -pe '\''s/.pt_(.*)_$(env).o/$(env)\/pt_$$1.o $@/g'\'' < .x.$@ > $@); \
	rm .x.$@; \
	exit $$rc'
endif
endif

ifdef kcc_deps
.d.%: .dkcc.%
	@$(FIXDEPS) -o $@ -i $< -e $(env) -v
endif

###############################################################################
#                        End of gts/src/config/makefile
###############################################################################
