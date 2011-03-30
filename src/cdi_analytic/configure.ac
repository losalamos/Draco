dnl-------------------------------------------------------------------------dnl
dnl draco/src/cdi_analytic configure.in
dnl Thomas M. Evans
dnl Fri Aug 24 12:28:30 2001
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for cdi_analytic
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([Analytic_Gray_Opacity.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(cdi_analytic/config.h:config.h.in)
package="cdi_analytic"

dnl VENDOR SETUPS
dnl AC_MPI_SETUP(pkg)
dnl AC_SHMEM_SETUP(pkg)
dnl AC_PCGLIB_SETUP(pkg)

AC_INSTALL_HEADERS
AC_INSTALL_LIB
AC_NEEDS_LIBS(cdi ds++)

AC_RUNTESTS(tstAnalytic_Gray_Opacity \
            tstnGray_Analytic_MultigroupOpacity tstAnalytic_EoS \
            tstAnalytic_Odfmg_Opacity, scalar)

dnl DRACO SETUP
AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-6_1_0)

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/cdi_analytic configure.in
dnl-------------------------------------------------------------------------dnl
