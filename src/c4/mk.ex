# Example C4 makefile (for dumb vendor make's).

# Set your compilation macros here.  Don't forget to add a -Iwherever
# to CXXFLAGS to point to your DS++ installation.  Also, note that
# config.h cues off CPP macros to determine the supporting message
# passing facility.  Your compiler driver may take care of this for
# you, but if not, add a -Dwhatever to CXXFLAGS.  Keep the -I.., since
# within C4, files are located via #include "c4/xxx.h", etc.

CXX = CC
CXXFLAGS = -I..
RANLIB	= ranlib

SRCS =	\
        NodeInfo.cc \
        SpinLock.cc \
        Sync.cc \
        global.cc \
        C4_Req.cc \
        BSwap.cc \
        DeepObj.cc \
        DeepBSwap.cc \
        ArraySwap.cc \
        DynArraySwap.cc


OBJS = $(SRCS:.cc=.o)

.SUFFIXES:
.SUFFIXES: .cc .o

.cc.o:
	$(CXX) $(CXXFLAGS) -c $*.cc

# Build the library

libc4.a : $(OBJS)
	-rm $@
	ar r $@ $(OBJS)
	$(RANLIB) $@

# Build the test programs.
