###############################################################################
# gts/src/ds++/makefile
# Geoffrey Furnish
# 17 January 1994
###############################################################################
# @> Rules for constructing the GTS data structures library, libds.a.
###############################################################################

pkg := ds++

# Go get the standard includeable makefile for GTS components.

include ../config/makefile

# Non template files, needed by all compilation environments.

#libsrcs	= \
#	charTable.cc \
#	String.cc
libsrcs	= Assert.cc String.cc

# Files which implement template classes, needed for C-front style
# compilers. 

tmplsrcs = \
	DynArray.cc \
	SP.cc \
	Sort.cc

# Get the library built.

include ../config/make.lib

# Now concern ourselves with the business of test programs.

notstlib := 1
tstsrcs	:= $(wildcard tst*.cc)

tstsrcs := $(tstsrcs:tstSort.cc=)

include ../config/make.tests

# Testing out some dejagnu support

include ../config/make.deja

# And now documentation.

docsrcs = Array.h Array.cc DynArray.h DynArray.cc \
	List.h Map.h Map.cc Stack.h Stack.cc String.h String.cc

include ../config/make.doc

# Go pick up the standard signature for makefiles.

cc_files        := $(libsrcs) $(tstsrcs)

include ../config/make.end

###############################################################################
# 			End of gts/src/ds++/makefile
###############################################################################
