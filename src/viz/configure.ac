dnl-------------------------------------------------------------------------dnl
dnl draco/src/viz configure.in
dnl Thomas M. Evans
dnl Fri Jan 21 16:29:46 2000
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for viz
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([Ensight_Translator.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(viz/config.h:config.h.in)
package="viz"

AC_INSTALL_HEADERS
AC_INSTALL_LIB
AC_NEEDS_LIBS(traits ds++)

dnl tests
AC_RUNTESTS(tstEnsight_Translator \
            tstEnsight_Stream, scalar)

dnl DRACO SETUP
AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-6_0_0)

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in \
                 test/cell_data:test/cell_data.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/viz configure.in
dnl-------------------------------------------------------------------------dnl
