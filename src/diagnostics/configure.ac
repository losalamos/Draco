dnl-------------------------------------------------------------------------dnl
dnl diagnostics/configure.ac
dnl Thomas M. Evans
dnl Fri Dec  9 11:05:13 2005
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl Copyright © 2003 The Regents of the University of California.
dnl-------------------------------------------------------------------------dnl
dnl configure.ac for diagnostics
dnl-------------------------------------------------------------------------dnl

AC_INIT

dnl As an argument to this macro provide the name of a file that
dnl exists in the local package directory.  
AC_CONFIG_SRCDIR([Diagnostics.hh])  

AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(diagnostics/config.h:config.h.in)
package="diagnostics"

dnl Draco version
AC_DRACO_PREREQ

AC_DEFINE(AUTOCONF autoconf)

dnl ----------------------------------------
dnl Vendor setups: 
dnl ----------------------------------------

AC_MPI_SETUP(pkg)

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

AC_NEEDS_LIBS(c4 ds++)

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

dnl AC_RUNTESTS(tstutils, 1 2)
AC_RUNTESTS(tstDiagnostics \
            tstTiming, 1)

dnl ----------------------------------------
dnl SETUP TIMING AND DIAGNOSTICS

AC_ARG_WITH(draco-diagnostics,
   [  --with-draco-diagnostics[=0-7]
                             set draco diagnostics level (default 0)])

AC_ARG_WITH(draco-timing,
   [  --with-draco-timing[=0-2]
                             set draco timing level (default 0)])

AC_MSG_CHECKING("for diagnostics")
# setup diagnostics level
if test "${with_draco_diagnostics:=default}" = no ; then
   with_draco_diagnostics='0'
elif test "${with_draco_diagnostics}" = yes ; then
   with_draco_diagnostics='7'
fi

if test "${with_draco_diagnostics}" != default ; then
   AC_DEFINE_UNQUOTED(DRACO_DIAGNOSTICS, $with_draco_diagnostics)
fi
AC_MSG_RESULT("$with_draco_diagnostics")

AC_MSG_CHECKING("for timing diagnostics")
# setup timing level
if test "${with_draco_timing:=default}" = no ; then
   with_draco_timing='0'
elif test "${with_draco_timing}" = yes ; then
   with_draco_timing='2'
fi

if test "${with_draco_timing}" != default ; then
   AC_DEFINE_UNQUOTED(DRACO_TIMING, $with_draco_timing)
fi
AC_MSG_RESULT("$with_draco_timing")

dnl ----------------------------------------

dnl ----------------------------------------
dnl Finish Draco setup:
dnl ----------------------------------------

AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-6_0_0)

dnl Create Makefiles from config templates.

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in])

dnl Create test/Makefile.options, if it exists.

AC_CHECK_FILE($srcdir/test/Makefile.options.in, \
AC_CONFIG_FILES([test/Makefile.options:test/Makefile.options.in]))

AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of diagnostics/configure.ac
dnl-------------------------------------------------------------------------dnl
