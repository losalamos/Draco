###############################################################################
# gts/src/nml/makefile
# Geoffrey Furnish
# 25 January 1994
###############################################################################
# @> Rules for building the Physical Dynamics Namelist Utility.
###############################################################################

pkg := nml

# Go get the standard includeable makefile for GTS components.

need_tcltk	= yes

include ../config/makefile

# Files needed in all compilation environments.

LIBSRCS	= \
	Item.cc \
	Block.cc \
	Group.cc \
	Items.cc \
	nmostream.cc \
	nmistream.cc

# Explicit template instantiation files.

ifdef preinstantiate

pt_files = $(wildcard *.pt)
ptcc_files = $(pt_files:%.pt=.pt_%_$(env).cc)

endif

TKLIBSRCS = $(LIBSRCS:Group.cc=tkGroup.cc) \
	nmitcl.cc \
	TraceInfo.cc

#ifndef tk4
#TKLIBSRCS += tkMain.cc
#endif

SRCS	= $(LIBSRCS) tkMain.cc nmitcl.cc TraceInfo.cc nmtest.cc tkAppInit.cc

ifdef parallel
spartan = yes
endif

ifdef spartan
ALL_SRCS	= $(LIBSRCS)
else
# tkGroup.cc is fake, this case handled by an explicit auto dependency
# rule. See below.

ALL_SRCS	= $(SRCS) tkGroup.cc
endif

ifdef use_gnu_ar_rule
NMLLIB_OBJS	= $(LIBSRCS:%.cc=libnml_$(env).a(%.o))
TKNMLLIB_OBJS	= $(TKLIBSRCS:%.cc=libnmltk_$(env).a(%.o))
else
NMLLIB_OBJS	= $(LIBSRCS:%.cc=$(env)/%.o)
TKNMLLIB_OBJS	= $(TKLIBSRCS:%.cc=$(env)/%.o)
endif

ifdef preinstantiate
NMLLIB_OBJS	+= $(ptcc_files:.%_$(env).cc=$(env)/%.o)
TKNMLLIB_OBJS	+= $(ptcc_files:.%_$(env).cc=$(env)/%.o)
endif

libnml		= libnml_$(env).$(lib_suffix)
libnmltk	= libnmltk_$(env).$(lib_suffix)

ifdef spartan
default: deps $(libnml)
else
default: deps $(libnml) $(libnmltk) index
tests: default nmtest_$(env) nml_wish_$(env)
endif

deps:
ifneq ($(dodeps),no)
	(cd ../ds++;   $(MAKE) $(subnj))
endif

resolve_libs	:= -L../ds++ -lds++_$(env)

#-----------------------------------------------------------------------------#
# There is a uniquely infuriating problem with HP/UX, whereby putting
# a new .o into a large archive library takes forever and a day.
# Although I really like the gmake archive rule, I'm gonna toss it on
# HP where ar sucks so bad.

ifdef use_gnu_ar_rule
libnml_$(env).a: $(NMLLIB_OBJS)
	$(ranlib) $@
else
libnml_$(env).a : $(NMLLIB_OBJS)
	-rm -f $@
	$(lib_cmd) $@ $(NMLLIB_OBJS) $(resolve_libs)

#	$(ar) qc $@ $(NMLLIB_OBJS)
#	$(ranlib) $@
endif

libnml_$(env).$(shlib_suffix) : $(NMLLIB_OBJS)
	$(shlib_cmd) $@ $(NMLLIB_OBJS)

ifdef use_gnu_ar_rule
libnmltk_$(env).a: $(TKNMLLIB_OBJS)
	$(ranlib) $@
else
# Introduce the dependence of libnmltk upon libnml in order to prevent
# their template closure phases from running concurrently if building
# on an SMP.

libnmltk_$(env).a : libnml_$(env).a $(TKNMLLIB_OBJS)
	-rm -f $@
	$(lib_cmd) $@ $(TKNMLLIB_OBJS)

#	$(ar) qc $@ $(TKNMLLIB_OBJS)
#	$(ranlib) $@
endif

libnmltk_$(env).$(shlib_suffix) : $(TKNMLLIB_OBJS)
	$(shlib_cmd) $@ $(TKNMLLIB_OBJS)

$(env)/tkGroup.o : Group.cc
	$(C++) $(C++FLAGS) -DTK -c Group.cc -o $@ 

nml_library = $(shell cd ..; pwd)/nml

$(env)/nmitcl.o : nmitcl.cc
	$(CXX) -c $(CXXFLAGS) -DNML_LIBRARY=\"$(nml_library)\" $< -o $@

nmtest_$(env): $(env)/nmtest.o $(dep_libnml) $(dep_libds++)
	$(C++) $(C++FLAGS) -o $@ $(env)/nmtest.o -L. -lnml_$(env) \
		$(LINK_LIBS)

nml_wish_$(env): $(env)/tkAppInit.o $(dep_libnmltk) $(dep_libds++)
	$(C++) $(C++FLAGS) -o $@ $(env)/tkAppInit.o -L. -lnmltk_$(env) \
		$(LINK_LIBS)

index: tclIndex

tcl_files = $(wildcard *.tcl *.itk *.ith *.itm )

tclIndex: $(tcl_files)
	$(INDEX)

.d.tkGroup_$(env): Group.cc
	@$(SHELL) -ec '$(MAKEDEPEND) \
	$(CXX_DEPEND_FLAG) -DTK $< | \
	sed '\''s/Group.o/tk& $@/g'\'' > $@'

# Go pick up the standard signature for makefiles.

cc_files        := $(ALL_SRCS)

nml_xcld        := init.nml

include ../config/make.end

# End of gts/src/nml/makefile
