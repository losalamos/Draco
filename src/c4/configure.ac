dnl-------------------------------------------------------------------------dnl
dnl draco/src/c4 configure.in
dnl Thomas M. Evans
dnl Fri Jan  8 14:46:27 1999
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for ds++
dnl-------------------------------------------------------------------------dnl
dnl $Id$
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([C4_Req.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(c4/config.h:config.h.in)
package="c4"

dnl VENDOR SETUPS
AC_MPI_SETUP(pkg)

AC_INSTALL_HEADERS
AC_INSTALL_LIB
AC_NEEDS_LIBS(ds++)

dnl DRACO SETUP

AC_RUNTESTS( \
            tstC4_Req \
            tstTime \
            tstReduction \
            tstPingPong \
            tstComm_Dup \
            tstBroadcast \
            tstAbort \
            tstScalar \
            , 1 2 4)
AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-5_11_0)

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/c4 configure.in
dnl-------------------------------------------------------------------------dnl
