dnl-------------------------------------------------------------------------dnl
dnl draco/src/xm configure.in
dnl Randy M. Roberts
dnl Mon Jan 11 10:07:10 MST 1999
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for xm
dnl-------------------------------------------------------------------------dnl
dnl $Id$
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([MT_traits.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(traits/config.h:config.h.in)
package="traits"

AC_INSTALL_HEADERS
AC_INSTALL_LIB

dnl draco libraries
AC_NEEDS_LIBS(ds++)

dnl runtesting
AC_RUNTESTS(tstViz_Traits, scalar)

AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-6_0_0)

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/c4 configure.in
dnl-------------------------------------------------------------------------dnl
