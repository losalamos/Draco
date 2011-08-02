dnl-------------------------------------------------------------------------dnl
dnl units/configure.ac
dnl Kelly Thompson
dnl Mon Oct 27 10:20:33 2003
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for units
dnl-------------------------------------------------------------------------dnl
dnl $Id$
dnl-------------------------------------------------------------------------dnl

AC_INIT
dnl AC_CONFIG_SRCDIR([Release.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(units/config.h:config.h.in)
package="units"

AC_INSTALL_HEADERS
AC_INSTALL_LIB
AC_NEEDS_LIBS(ds++)

AC_RUNTESTS(tstUnitSystemEnums tstFundUnit tstUnitSystemType tstUnitSystem tstPhysicalConstants, scalar)

AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-6_2_0)

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in
                 test/Makefile:../../config/Makefile.test.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of units/configure.ac
dnl-------------------------------------------------------------------------dnl
