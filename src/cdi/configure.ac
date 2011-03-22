dnl-------------------------------------------------------------------------dnl
dnl draco/src/cdi configure.in
dnl Kelly Thompson
dnl Thu Jun 22 16:18:14 2000
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for cdi
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([CDI.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(cdi/config.h:config.h.in)
package="cdi"

AC_INSTALL_HEADERS
AC_INSTALL_LIB
AC_NEEDS_LIBS(ds++)

AC_RUNTESTS(tDummyOpacity tDummyEoS tCDI, scalar) 

dnl DRACO SETUP
AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-6_1_0)

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/cdi configure.in
dnl-------------------------------------------------------------------------dnl
