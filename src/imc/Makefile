###############################################################################
# imctest makefile
# Thomas M. Evans
###############################################################################
# @> IMCTEST package makefile
###############################################################################

##---------------------------------------------------------------------------##
## SETTINGS
##---------------------------------------------------------------------------##

CC = KCC
INCLUDE = -I../ds++ 
LIB = -L../ds++ -lds++_sunkaig -lm
OPTIONS = --strict -x 
COPTIONS = -M
RMREC = rm -r
RM = rm

##---------------------------------------------------------------------------##
## OBJECTS
##---------------------------------------------------------------------------##

objects = Layout.o OS_Builder.o OS_Mesh.o OS_Parser.o Opacity.o \
          Particle.o Random.o XYCoord_sys.o XYZCoord_sys.o \
          Opacity_Builder.o Mat_State.o Opacity.o

sources = Layout.cc OS_Builder.cc OS_Mesh.cc OS_Parser.cc Opacity.cc \
          Particle.cc Random.cc XYCoord_sys.cc XYZCoord_sys.cc \
          Opacity_Builder.cc Mat_State.cc Opacity.cc

##---------------------------------------------------------------------------##
## OBJECT FILE INCLUDE DEPENDENCIES
##---------------------------------------------------------------------------##

%.cd: %.cc %.hh
	$(CC) $< $(INCLUDE) -M > $@.t
	nawk '/.hh/ {print $2}' $@.t > $@
	$(RM) $@.t

include $(sources:.cc=.cd)

##---------------------------------------------------------------------------##
## COMPILE TARGETS FOR CLASS OBJECTS
##---------------------------------------------------------------------------##

all: $(objects)

%.o: %.cc
	$(CC) -c $< $(INCLUDE) $(OPTIONS)

##---------------------------------------------------------------------------##
## CLEAN UTILITIES
##---------------------------------------------------------------------------##

cleanall: cleanobj cleankcc cleantest

cleanobj:
	$(RM) *.o

cleankcc:
	$(RM) *.ii

cleantest:
	$(RM) test*

cleandep:
	$(RM) *.cd

##---------------------------------------------------------------------------##
## TEST DRIVERS
##---------------------------------------------------------------------------##

# find test executable drivers
vpath test%.cc test

testbobj = testb.o OS_Builder.o OS_Mesh.o XYCoord_sys.o XYZCoord_sys.o \
           Layout.o Random.o OS_Parser.o Opacity_Builder.o Mat_State.o \
           Opacity.o

testbdep = OS_Parser.hh OS_Builder.hh OS_Mesh.hh Opacity_Builder.hh \
           Mat_State.hh Opacity.hh

testb.o: $(testbdep)

testb: $(testbobj)
	$(CC) -o $@ $(testbobj) $(LIB) $(OPTIONS)

###############################################################################
#			   end of imctest makefile
###############################################################################
