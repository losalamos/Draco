dnl-------------------------------------------------------------------------dnl
dnl shared_lib/configure.ac
dnl Thomas M. Evans
dnl Wed Apr 21 14:16:24 2004
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl Copyright © 2003 The Regents of the University of California.
dnl-------------------------------------------------------------------------dnl
dnl configure.in for shared_lib
dnl-------------------------------------------------------------------------dnl

AC_INIT

dnl As an argument to this macro provide the name of a file that
dnl exists in the local package directory.  
AC_CONFIG_SRCDIR([Shared_Lib.hh])  

AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(shared_lib/config.h:config.h.in)
package="shared_lib"

dnl Draco version
AC_DRACO_PREREQ

dnl ----------------------------------------
dnl Vendor setups: 
dnl ----------------------------------------

dnl Uncomment any vendor setups that are required for this package.
dnl For example, if this package has parallel programming directives
dnl you need to uncomment AC_MPI_SETUP.

dnl AC_MPI_SETUP(pkg)
dnl AC_SHMEM_SETUP(pkg)
dnl AC_SPRNG_SETUP(pkg)
dnl AC_PCGLIB_SETUP(pkg)
dnl AC_TRILINOS_SETUP(pkg)
dnl AC_LAPACK_SETUP(pkg)
dnl AC_PARMETIS_SETUP(pkg)
AC_DLOPEN_SETUP(pkg)

dnl This macro should be uncommented if the current package is a Draco
dnl client (i.e.: This package exists in Capsaicin, Milagro, etc.)

dnl AC_DRACO_VENDOR_SETUP

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

AC_NEEDS_LIBS(ds++)

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

AC_RUNTESTS(tstShared_Lib, scalar)

dnl ----------------------------------------
dnl Finish Draco setup:
dnl ----------------------------------------

AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-6_1_0)

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of shared_lib/configure.ac
dnl-------------------------------------------------------------------------dnl
