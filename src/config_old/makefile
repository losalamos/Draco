###############################################################################
# gts/src/config/makefile
# Geoffrey Furnish
# 26 June 1993
###############################################################################
# @> This "makefile" is included by the makefiles in each component
# @> directory.  The intention is that this makefile will hold all the
# @> definitions required to configure compilers and tools for building
# @> the GTS components.  The makefile for each component will only need
# @> to include this one and then do its own component specific stuff
# @> without having to worry about how to set up all the compile and link
# @> options. 
###############################################################################

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

# Okay, we're gonna start a new thing here...  All basic libs should
# be under TOOL_DIR by default, but we'll still let the guy override
# if needed.

TOOL_DIR= $(HOME)/gts

GTS_SUPPORT_ENV_LIBDIR	= $(TOOL_DIR)/lib
GTS_SUPPORT_ENV_INCDIR	= $(TOOL_DIR)/include
GTS_MAN_DIR		= $(TOOL_DIR)/man

ITSH  = $(TOOL_DIR)/bin/itcl_sh
MKINDEX	= $(TOOL_DIR)/library/itcl_mkindex.tcl
INDEX = echo 'source $(MKINDEX) ; itcl_mkindex . *.tcl ; exit' | $(ITSH)

itk	= yes	# Hardwire it on for a while, see how we like the new stuff.

ifdef itk
ITCLSH = $(TOOL_DIR)/bin/itclsh
INDEX = echo 'auto_mkindex . *.tcl *.itk *.ith *.itm ; exit' | $(ITCLSH)
endif

# Set up list of all valid environmentals

envs = tgcc elf mpi hp pgn

# Ask not to inherit certain variables.

unexport need_gtslib need_gf++ need_tcltk need_namelist need_plplot need_x11 need_python

#-----------------------------------------------------------------------------#
# Configure for TCL/Tk

TK_INCS		= -I$(TOOL_DIR)/include
TK_LIBS		= -L$(TOOL_DIR)/lib -litk -litcl -lBLT -ldpnetwork -ltk -ltcl

# Try it this way, GTS_SUPPORT_ENV_LIBDIR assumed.
TK_LIBS		= -litk -litcl -lBLT -ldpnetwork -ltk -ltcl
TK_LIBS		= -litk -litcl -lBLT -ltk -ltcl

#-----------------------------------------------------------------------------#
# Configure for Physical Dynamics C++ Namelist Utility.

NML_LIB_DIR	= $(GTS_SRC_TOP)/nml

ifeq ($(need_tcltk),yes)
NML_LIBS	= -L$(NML_LIB_DIR) -lnmltk_$(env)
else
NML_LIBS	= -L$(NML_LIB_DIR) -lnml_$(env)
endif

#-----------------------------------------------------------------------------#
# Configure for PLPLOT.

plplotf		= -lplplotf
plplotftk	= -lplplotftk

ifeq ($(need_tcltk),yes)
PLPLOT_LIBS	= -L$(PLPLOT_LIB_DIR) $(plplotftk)
else
PLPLOT_LIBS	= -L$(PLPLOT_LIB_DIR) $(plplotf)
endif

# Try it this way, GTS_SUPPORT_ENV_LIBDIR assumed.

ifeq ($(need_tcltk),yes)
PLPLOT_LIBS	= $(plplotftk)
else
PLPLOT_LIBS	= $(plplotf)
endif

ifeq (${HOSTTYPE},hp9000s700)	# HP/UX
FFT_LIB		= -L$(HOME)/lib -lfftf -lfutilsf
ifeq ($(OC),on)
FFT_LIB		= -L$(HOME)/lib -lfftf -lfutilsf -lcl
endif
ifeq ($(cxx),gcc)
FFT_LIB		= -L$(HOME)/lib -lfftf -lfutilsf -lcl
endif
endif

ifeq (${HOSTTYPE},i386)		# Linux
TK_LIBS		= -L$(TOOL_DIR)/lib -litcl -lBLT -ltk -ltcl
TK_LIBS		= -litk -litcl -lBLT -ltk -ltcl

# Planning to keep utils stuff in ~/lib, but most of it doesn't even
# work with f2c on Linux.  Right now only libfft.a is checked out.

FFT_LIB		= -L$(HOME)/lib -lfftf -lfutilsf -lf2c
FFT_LIB		= -L$(TOOL_DIR)/lib -lfftf -lfutilsf -lf2c
FFT_LIB		= -L$(GTS_SRC_TOP)/fft -lfftf_$(env)
endif

#-----------------------------------------------------------------------------#
# Configure for X11

X11_LIBS	= -lX11 -lm

ifeq (${HOSTTYPE},sparc)	# MCC
TK_INCS		+= -I/usr/local/X11R5/include
X11_LIBS	= -L/usr/local/X11R5/lib -lX11 -lm
FFT_LIB		= -L$(HOME)/lib -lfftf -lfutilsf -lF77 -lM77
endif

ifeq (${HOSTTYPE},i386)		# Linux
X11_LIBS	= -L/usr/X11R6/lib -lX11 -lm
endif

#-----------------------------------------------------------------------------#
# Configure for debugging malloc

ifeq ($(dbmalloc),on)
DBMALLOC_LIB = -L$(TOOL_DIR)/lib -lmalloc
DBMALLOC_LIB = -lmalloc
else
DBMALLOC_LIB =
endif

#-----------------------------------------------------------------------------#
# Configure for C4 (Cannonical Classes for Concurrency Control)

C4_LIB		= -L$(GTS_SRC_TOP)/c4 -lc4_$(env)

#-----------------------------------------------------------------------------#
# Configure for GTS development.

ifndef GTS_SRC_TOP
GTS_SRC_TOP	= ..
endif

TRACER_LIB	= -L$(GTS_SRC_TOP)/tracer -ltracer_$(env)
DS++_LIB	= -L$(GTS_SRC_TOP)/ds++ -lds++_$(env)
GTS_LIB		= -L$(GTS_SRC_TOP)/base -lgts_$(env)
PTL_LIB		= -L$(GTS_SRC_TOP)/ptl -lptl_$(env)
GF++_LIBS	= -L$(GTS_SRC_TOP)/gf++ -lgf++_$(env)
SPDF_LIB	= -L$(GTS_SRC_TOP)/spdf -lspdf_$(env)
UTIL_LIB	= -L$(GTS_SRC_TOP)/util -lutil_$(env)
LINALG_LIB	= -L$(GTS_SRC_TOP)/linalg -llinalg_$(env)
FLDEQNS_LIB	= -L$(GTS_SRC_TOP)/fldeqns -lfldeqns_$(env)

#-----------------------------------------------------------------------------#
# Set up some tools.

MAKEDEPEND	:= gcc -MM -U__GNUC__

NMLGEN		:= $(GTS_SRC_TOP)/../tools/nmlgen

PERSTALYZER	:= $(GTS_SRC_TOP)/../tools/perstalyzer

INSTANTIALYZER	:= $(GTS_SRC_TOP)/../tools/instantialyzer

AUTODOC		:= $(GTS_SRC_TOP)/../tools/autodoc

FIXDEPS		:= $(GTS_SRC_TOP)/../tools/fixdeps

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

# Cut out some makefile clutter.

ifeq ($(PLPLOT_INCS),$(TK_INCS))
TK_INCS =
endif

ar	= ar
ranlib	= ranlib

# Set this to -DNO_XCPT if Assert.h can't detect this compiler for
# some reason, or other inanity occurs.

XCPT_FLAG =

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

oscfgfile = $(GTS_SRC_TOP)/config/$(uname).cf

ifeq ($(wildcard $(oscfgfile)),$(oscfgfile))
include $(oscfgfile)
endif

hstcfgfile = $(GTS_SRC_TOP)/config/$(host).cf

ifeq ($(wildcard $(hstcfgfile)),$(hstcfgfile))
include $(hstcfgfile)
endif

###############################################################################

# junk.1 was here.

###############################################################################
# Now setup the main compilation control variables.

# For now, just put all conceivable -I's and so forth on the compile
# line so that we can hopefully minimize the Cfront template
# instantiation trauma.

# This used to be in CXX_INC_FLAGS: -I$(GTS_SRC_TOP)/ds++ 

CXX_INC_FLAGS = -I$(GTS_SRC_TOP) \
	-I$(GTS_SUPPORT_ENV_INCDIR) $(PY_INC_FLAG)

CXXFLAGS = $(SYS_FLAGS) $(CXX_INC_FLAGS) $(DEBUG_FLAG) $(XCPT_FLAG) $(ASSERT_FLAG) $(PYTHON_FLAG)

# Need these because I have some cases of explicit rules in my various
# makefiles, written back when I was using gmake before the
# conversion from C++ to CXX, etc.

C++ = $(CXX)
C++FLAGS = $(CXXFLAGS)

#-----------------------------------------------------------------------------#
# This section is coded in this manner b/c for example, the user might
# (conceivably) want Tcl/Tk without PLplot.  In other words, we don't
# fold the logic for need_tcltk directly into that for plplot, to
# allow more flexibility of choice to the user.

ifeq ($(need_tracer),yes)
need_complex	= yes
need_gf++	= yes
endif

ifeq ($(need_gf++),yes)
need_gtslib	= yes
need_plplot	= yes
endif

ifeq ($(need_ptl),yes)
need_gtslib	= yes
endif

ifeq ($(need_gtslib),yes)
need_namelist	= yes
need_c4		= yes
need_spdf	= yes
need_util	= yes
need_fft	= yes
endif

ifeq ($(need_plplot),yes)
need_tcltk	= yes
need_tooldir	= yes
endif

ifeq ($(need_tcltk),yes)
need_x11	= yes
need_tooldir	= yes
endif

#-----------------------------------------------------------------------------#
# Now construct the macro used for linking programs.

LINK_LIBS	:=
ifdef SYS_LINK_FLAGS
LINK_LIBS	+= $(SYS_LINK_FLAGS)
endif
ifeq ($(need_tracer),yes)
LINK_LIBS	+= $(TRACER_LIB)
endif
ifeq ($(need_gf++),yes)
LINK_LIBS	+= $(GF++_LIBS)
endif
ifeq ($(need_fldeqns),yes)
LINK_LIBS	+= $(FLDEQNS_LIB)
endif
ifeq ($(need_ptl),yes)
LINK_LIBS	+= $(PTL_LIB)
endif
ifeq ($(need_gtslib),yes)
LINK_LIBS	+= $(GTS_LIB)
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
ifeq ($(need_linalg),yes)
LINK_LIBS	+= $(LINALG_LIB)
endif
ifeq ($(need_c4),yes)
LINK_LIBS	+= $(C4_LIB)
endif
LINK_LIBS	+= $(DS++_LIB)
ifeq ($(need_tooldir),yes)
LINK_LIBS	+= -L$(GTS_SUPPORT_ENV_LIBDIR)
endif
ifeq ($(need_python),yes)
LINK_LIBS	+= $(PYLIBS)
endif
ifeq ($(need_plplot),yes)
LINK_LIBS	+= $(PLPLOT_LIBS)
endif
ifeq ($(need_tcltk),yes)
LINK_LIBS	+= $(TK_LIBS)
endif
ifeq ($(need_x11),yes)
LINK_LIBS	+= $(X11_LIBS)
endif
ifeq ($(need_fft),yes)
LINK_LIBS	+= $(FFT_LIB)
endif
# dbmalloc support needs to be added
LINK_LIBS	+= $(DBMALLOC_LIB)
ifdef MPI_LIB
LINK_LIBS	+= $(MPI_LIB)
endif
ifdef GCC_LIBS
LINK_LIBS	+= $(GCC_LIBS)
endif
ifeq ($(need_complex),yes)
LINK_LIBS	+= $(COMPLEX_LIB)
endif
ifdef SYS_LIBS
LINK_LIBS	+= $(SYS_LIBS)
endif
LINK_LIBS	+= -lm

#-----------------------------------------------------------------------------#
# Now construct some macros for use in writing dependencies upon libraries.

ifdef make_shared_lib
lib_suffix	= $(shlib_suffix)
else
lib_suffix	= a
endif

dep_libtracer	= $(GTS_SRC_TOP)/tracer/libtracer_$(env).$(lib_suffix)
dep_libds++	= $(GTS_SRC_TOP)/ds++/libds++_$(env).$(lib_suffix)
dep_libnml	= $(GTS_SRC_TOP)/nml/libnml_$(env).$(lib_suffix)
dep_libnmltk	= $(GTS_SRC_TOP)/nml/libnmltk_$(env).$(lib_suffix)
dep_libc4	= $(GTS_SRC_TOP)/c4/libc4_$(env).$(lib_suffix)
dep_libspdf	= $(GTS_SRC_TOP)/spdf/libspdf_$(env).$(lib_suffix)
dep_libutil	= $(GTS_SRC_TOP)/util/libutil_$(env).$(lib_suffix)
dep_libptl	= $(GTS_SRC_TOP)/ptl/libptl_$(env).$(lib_suffix)
dep_libgts	= $(GTS_SRC_TOP)/base/libgts_$(env).$(lib_suffix)
dep_libgf++	= $(GTS_SRC_TOP)/gf++/libgf++_$(env).$(lib_suffix)
dep_libgtp	= $(GTS_SRC_TOP)/gtp/libgtp_$(env).$(lib_suffix)

dep_gtsmod	= $(GTS_SRC_TOP)/pymod/libgtsmodule_$(env).a

#-----------------------------------------------------------------------------#
# Now finally, some implicit rules which we will be using constantly,
# which depend critically on the value of $(env).

$(env)/%.o : %.cc
	$(CXX) -c $(CXXFLAGS) $*.cc -o $@

.pt_%_$(env).cc : %.pt
	$(INSTANTIALYZER) -o $@ $*.pt

$(env)/%.o : .%_$(env).cc
	$(CXX) -c $(CXXFLAGS) $(usrflg) .$*_$(env).cc -o $@

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

ifdef kcc_deps
.d.%: .dkcc.%
	@$(FIXDEPS) -o $@ -i $< -e $(env) -v
endif

###############################################################################
#                        End of gts/src/config/makefile
###############################################################################
