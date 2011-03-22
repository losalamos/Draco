dnl-------------------------------------------------------------------------dnl
dnl draco/src/lapack_wrap configure.in
dnl Thomas M. Evans
dnl Thu Aug 29 11:06:46 2002
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for lapack_wrap
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([Blas.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(lapack_wrap/config.h:config.h.in)
package="lapack_wrap"

dnl VENDOR SETUPS
AC_LAPACK_SETUP(pkg)

AC_INSTALL_HEADERS
AC_INSTALL_LIB
AC_NEEDS_LIBS(ds++)

# DEFINE BLAS/LAPACK FUNCTIONS
AC_DEFINE(LW_AUTOCONF)
AC_DEFINE(SCOPY, scopy_)
AC_DEFINE(DCOPY, dcopy_)
AC_DEFINE(SSCAL, sscal_)
AC_DEFINE(DSCAL, dscal_)
AC_DEFINE(SDOT, sdot_)
AC_DEFINE(DDOT, ddot_)
AC_DEFINE(DAXPY, daxpy_)
AC_DEFINE(SAXPY, saxpy_)
AC_DEFINE(SNRM2, snrm2_)
AC_DEFINE(DNRM2, dnrm2_)

# testing setup
AC_RUNTESTS(tstBlas_Level_1, scalar)

dnl DRACO SETUP
AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-6_1_0)

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/lapack_wrap configure.in
dnl-------------------------------------------------------------------------dnl
