# SETTINGS

CC = KCC
INCLUDE = -I../ds++ 
LIB = -L../ds++ -lds++_sunkaig -lm
OPTIONS = --strict -x 
COPTIONS = -M
RM = rm -r

# MAIN DRIVER PROGRAM

all: Particle.o XYCoord_sys.o XYZCoord_sys.o Layout.o OS_Mesh.o \
     Random.o OS_Builder.o Random.o Opacity.o OS_Parser.o

# CLASS MODULES

Particle.o: Particle.cc Particle.hh Coord_sys.hh Global.hh OS_Mesh.hh \
            Names.hh Random.hh Opacity.hh
	${CC} -c Particle.cc ${INCLUDE} ${OPTIONS}

XYCoord_sys.o: XYCoord_sys.cc XYCoord_sys.hh Coord_sys.hh Names.hh \
               Global.hh
	${CC} -c XYCoord_sys.cc ${INCLUDE} ${OPTIONS}

XYZCoord_sys.o: XYZCoord_sys.cc XYZCoord_sys.hh Coord_sys.hh Names.hh \
                Global.hh
	${CC} -c XYZCoord_sys.cc ${INCLUDE} ${OPTIONS}

Layout.o: Layout.cc Layout.hh Names.hh
	${CC} -c Layout.cc ${INCLUDE} ${OPTIONS}

OS_Mesh.o: OS_Mesh.cc OS_Mesh.hh Layout.hh Coord_sys.hh Global.hh Names.hh
	${CC} -c OS_Mesh.cc ${INCLUDE} ${OPTIONS}

Random.o: Random.cc Random.hh Global.hh Names.hh
	${CC} -c Random.cc ${INCLUDE} ${OPTIONS}

Opacity.o: Opacity.cc Opacity.hh OS_Mesh.hh Layout.hh Names.hh
	${CC} -c Opacity.cc ${INCLUDE} ${OPTIONS}

OS_Builder.o: OS_Builder.cc OS_Builder.hh OS_Mesh.hh Layout.hh \
              Coord_sys.hh Names.hh XYCoord_sys.hh XYZCoord_sys.hh
	${CC} -c OS_Builder.cc ${INCLUDE} ${OPTIONS}

OS_Parser.o: OS_Parser.cc OS_Parser.hh OS_Mesh.hh Names.hh
	${CC} -c OS_Parser.cc ${INCLUDE} ${OPTIONS}

OS_Parser.o:

# UTILITIES

cleanall:
	${RM} *.o

cleankcc:
	${RM} *.ii

# TEST DRIVERS

testm: testm.o XYCoord_sys.o XYZCoord_sys.o OS_Mesh.o Layout.o Random.o
	${CC} -o testm testm.o XYCoord_sys.o XYZCoord_sys.o \
                 OS_Mesh.o Layout.o Random.o ${LIB} ${OPTIONS}

testm.o: testm.cc Global.hh Layout.hh XYCoord_sys.hh XYZCoord_sys.hh \
         Random.hh
	${CC} -c testm.cc ${INCLUDE} ${OPTIONS}

testb: testb.o OS_Builder.o OS_Mesh.o XYCoord_sys.o XYZCoord_sys.o \
       Layout.o Random.o OS_Parser.o
	${CC} -o testb testb.o OS_Builder.o OS_Mesh.o XYCoord_sys.o \
                 XYZCoord_sys.o Layout.o Random.o OS_Parser.o ${LIB} ${OPTIONS}

testb.o: testb.cc OS_Parser.hh OS_Builder.hh OS_Mesh.hh
	${CC} -c testb.cc ${INCLUDE} ${OPTIONS}
