dnl-------------------------------------------------------------------------dnl
dnl draco/src/timestep configure.in
dnl Randy M. Roberts
dnl 
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for timestep
dnl-------------------------------------------------------------------------dnl
dnl $Id$
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([ts_manager.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(timestep/config.h:config.h.in)
package="timestep"

dnl VENDOR SETUPS
AC_MPI_SETUP(pkg)

AC_INSTALL_HEADERS
AC_INSTALL_LIB
AC_NEEDS_LIBS(c4 ds++)

dnl DRACO SETUP

AC_RUNTESTS(tstTimeStep, 1 2)
AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-6_2_0)

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/timestep configure.in
dnl-------------------------------------------------------------------------dnl
