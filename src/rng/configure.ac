dnl-------------------------------------------------------------------------dnl
dnl draco/src/rng configure.in
dnl Thomas M. Evans
dnl Mon Jan 11 10:21:40 1999
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for rng
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([Release.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(rng/config.h:config.h.in)
package="rng"

# Toggle between autoconf and cmake build systems when constructing
# config.h from config.h.in
AC_DEFINE(AUTOCONF autoconf)

AC_LAPACK_SETUP(pkg)
AC_GSL_SETUP(pkg)

AC_INSTALL_HEADERS
AC_INSTALL_LIB
AC_NEEDS_LIBS(ds++)

AC_ARG_ENABLE([rng-nr],AS_HELP_STRING([--enable-rng-nr],
                                      [use non-reproducible random number streams]),
AC_DEFINE([RNG_NR]))   

# TESTS
AC_RUNTESTS(tstRnd_Control_Inline \
            tstSubrandom_Sequence \
            , scalar)

# DRACO SETUP
AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-5_21_0)

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/rng configure.in
dnl-------------------------------------------------------------------------dnl
