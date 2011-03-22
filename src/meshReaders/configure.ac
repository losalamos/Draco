dnl-------------------------------------------------------------------------dnl
dnl draco/src/meshReaders configure.in
dnl Todd Adams
dnl Tue Sep 21 15:09:59 1999
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for meshReaders
dnl-------------------------------------------------------------------------dnl

AC_INIT

dnl As an argument to this macro provide the name of a file that
dnl exists in the local package directory.  
AC_CONFIG_SRCDIR([Mesh_Reader.hh])

AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(meshReaders/config.h:config.h.in)
package="meshReaders"

dnl Draco version
AC_DRACO_PREREQ

dnl ----------------------------------------
dnl Install options:
dnl ----------------------------------------
AC_INSTALL_HEADERS
AC_INSTALL_LIB

dnl ----------------------------------------
dnl Package dependencies:
dnl ----------------------------------------

dnl List local package dependencies here as arguments to
dnl AC_NEEDS_LIBS.  If this is a Draco client (i.e.: a part of
dnl Capsaicin, Milagro, etc.), the dependencies on Draco libraries
dnl should be listed as arguments to AC_NEEDS_LIBS_DRACO.

AC_NEEDS_LIBS(mesh_element \
              ds++)

dnl ----------------------------------------
dnl Setup for unit tests:
dnl ----------------------------------------

dnl List local package dependencies for unit tests here as arguments
dnl to AC_NEEDS_LIBS_TEST.  Unit test dependencies are usually
dnl different than package dependencies.  If this is a Draco client,
dnl the unit test dependencies on Draco libraries should be listed as
dnl arguments to AC_NEEDS_LIBS_TEST_DRACO.

dnl AC_NEEDS_LIBS_TEST(c4 ds++)
dnl AC_NEEDS_LIBS_TEST_DRACO(c4 ds++)

dnl List the unit tests that should be run and the number of
dnl processors that should be used when running the test. If this is a
dnl parallel package the number of processors should be a an integer.
dnl If this is a scalar package, the number of processors required
dnl should be listed as "scalar".

AC_RUNTESTS(TestHexMeshReader, scalar)

dnl ----------------------------------------
dnl Finish Draco setup:
dnl ----------------------------------------

AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-6_1_0)

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in   \
                 test/Makefile:../../config/Makefile.test.in \
	         test/cube.mesh:test/cube.mesh.in \
	         test/quad.mesh:test/quad.mesh.in \
 	         test/slab.mesh:test/slab.mesh.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/meshReaders configure.in
dnl-------------------------------------------------------------------------dnl
