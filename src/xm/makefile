###############################################################################
# kull/src/xm/makefile
# Geoffrey Furnish
# 5 September 1996
###############################################################################

pkg := xm

need_libm	:= yes

# Go get the standard includeable makefile for Kull components.

include ../config/makefile

default: deps tests

# Now build the test programs

notstlib = yes
tstsrcs	= tstUV.cc

include ../config/make.tests

deps:
ifneq ($(dodeps),no)
	(cd ../ds++; $(MAKE) $(subnj))
endif

# Get DejaGnu support

include ../config/make.deja

# Go pick up the standard signature for makefiles.

cc_files	:= $(libsrcs) $(tstsrcs)

include ../config/make.end

# End of kull/src/xm/makefile
