# SETTINGS

CC = KCC
INCLUDE = -I../ds++ 
LIB = -L../ds++ -lds++_sunkaig -lm
OPTIONS = --strict -x 
COPTIONS = -M
RM = rm -r

# OBJECTS

objects = Layout.o OS_Builder.o OS_Mesh.o OS_Parser.o Opacity.o \
          Particle.o Random.o XYCoord_sys.o XYZCoord_sys.o

sources = Layout.cc OS_Builder.cc OS_Mesh.cc OS_Parser.cc Opacity.cc \
          Particle.cc Random.cc XYCoord_sys.cc XYZCoord_sys.cc

# INCLUDE DEPENDENCIES

%.cd: %.cc %.hh
	$(CC) $< $(INCLUDE) -M > $@.t
	nawk '/.hh/ {print $2}' $@.t > $@
	$(RM) $@.t

include $(sources:.cc=.cd)

# MAIN DRIVER PROGRAM

all: $(objects)

# CLASS OBJECTS

%.o: %.cc
	$(CC) -c $< $(INCLUDE) $(OPTIONS)

# UTILITIES

cleanall: cleanobj cleankcc cleandep

cleandep:
	$(RM) *.cd

cleanobj:
	$(RM) *.o

cleankcc:
	$(RM) *.ii

cleantest:
	rm test*

# TEST DRIVERS

vpath test%.cc test

testsources = testb.cc testm.cc

testmobj = testm.o XYCoord_sys.o XYZCoord_sys.o OS_Mesh.o Layout.o\
           Random.o

testbobj = testb.o OS_Builder.o OS_Mesh.o XYCoord_sys.o XYZCoord_sys.o \
           Layout.o Random.o OS_Parser.o

%.td: %.cc
	$(CC) $< $(INCLUDE) -M > $@.t
	nawk '/.hh/ {print $2}' $@.t > $@
	$(RM) $@.t

include $(testsources:.cc=.td)

testm: $(testmobj)
	$(CC) -o $@ $(testmobj) $(LIB) $(OPTIONS)

testb: $(testbobj)
	$(CC) -o $@ $(testbobj) $(LIB) $(OPTIONS)