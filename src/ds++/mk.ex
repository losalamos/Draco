# Example DS++ makefile (for dumb vendor make's).

CXX = CC
CXXFLAGS =
RANLIB	= ranlib

SRCS =	Array.cc DynArray.cc DynTable.cc List.cc Map.cc \
	SArray.cc SP.cc Stack.cc String.cc slist.cc

OBJS = $(SRCS:.cc=.o)

.SUFFIXES:
.SUFFIXES: .cc .o

.cc.o:
	$(CXX) $(CXXFLAGS) -c $*.cc

libds++.a : $(OBJS)
	-rm $@
	ar r $@ $(OBJS)
	$(RANLIB) $@