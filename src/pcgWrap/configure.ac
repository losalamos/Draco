dnl-------------------------------------------------------------------------dnl
dnl draco/src/pcgWrap configure.in
dnl Randy M. Roberts
dnl 
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for pcgWrap
dnl-------------------------------------------------------------------------dnl
dnl $Id$
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([PCG_Ctrl.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(pcgWrap/config.h:config.h.in)
package="pcgWrap"

AC_PCG_SETUP(pkg) dnl through pcgWrap
AC_MPI_SETUP(pkg)    dnl through c4
AC_LAPACK_SETUP(pkg) 

AC_INSTALL_HEADERS
AC_INSTALL_LIB
AC_NEEDS_LIBS(ds++)
AC_NEEDS_LIBS_TEST(c4)

AC_RUNTESTS(TstPCG, 1 2 5)

AC_DRACO_ENV

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/pcgWrap configure.in
dnl-------------------------------------------------------------------------dnl
