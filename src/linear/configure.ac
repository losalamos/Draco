dnl-------------------------------------------------------------------------dnl
dnl linear/configure.ac
dnl Kent G. Budge
dnl Wed Nov 10 08:02:06 2010
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl Copyright (C) 2006 Los Alamos National Security, LLC
dnl-------------------------------------------------------------------------dnl
dnl configure.ac for linear
dnl-------------------------------------------------------------------------dnl

AC_INIT

dnl As an argument to this macro provide the name of a file that
dnl exists in the local package directory.  
AC_CONFIG_SRCDIR([README.linear])  

AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(linear/config.h:config.h.in)
package="linear"

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
dnl AC_GSL_SETUP(pkg)
dnl AC_BLACS_SETUP(pkg)
dnl AC_SCALAPACK_SETUP(pkg)
dnl AC_HYPRE_SETUP(pkg)

dnl This macro should be uncommented if the current package is a Draco
dnl client (i.e.: This package exists in Capsaicin, Milagro, etc.)

dnl AC_DRACO_VENDOR_SETUP

dnl ----------------------------------------
dnl Install options:
dnl ----------------------------------------

AC_INSTALL_HEADERS
AC_INSTALL_LIB
dnl AC_INSTALL_EXECUTABLE

dnl ----------------------------------------
dnl Package dependencies:
dnl ----------------------------------------

dnl List local package dependencies here as arguments to
dnl AC_NEEDS_LIBS.  If this is a Draco client (i.e.: a part of
dnl Capsaicin, Milagro, etc.), the dependencies on Draco libraries
dnl should be listed as arguments to AC_NEEDS_LIBS_DRACO.

AC_NEEDS_LIBS(ds++)
dnl AC_NEEDS_LIBS_DRACO(ds++)

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

dnl AC_RUNTESTS(tstlinear, 1 2)
AC_RUNTESTS(tstrotate \
            tstfnorm \
            tstrsolv \
            tstsvdcmp \
            tstsvbksb \
            tstqrdcmp \
            tstqr_unpack \
            tstbtridag \
            tstgaussj \
            tstludcmp \
            tsttqli \
            tsttred2 \
            tstqrupdt \
           , scalar)

dnl ----------------------------------------
dnl Finish Draco setup:
dnl ----------------------------------------

AC_DRACO_ENV
AC_DRACO_AUTODOC

dnl Create Makefiles from config templates.

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in])

dnl Create test/Makefile.options, if it exists.

AC_CHECK_FILE($srcdir/test/Makefile.options.in, \
AC_CONFIG_FILES([test/Makefile.options:test/Makefile.options.in]))

AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of linear/configure.ac
dnl-------------------------------------------------------------------------dnl
