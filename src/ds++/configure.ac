dnl-------------------------------------------------------------------------dnl
dnl draco/src/ds++ configure.in
dnl Thomas M. Evans
dnl Thu Jan  7 16:58:31 1999
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for ds++
dnl-------------------------------------------------------------------------dnl
dnl $Id$
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([SP.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(ds++/config.h:config.h.in)
package="ds++"

AC_INSTALL_HEADERS
AC_INSTALL_LIB

AC_RUNTESTS(tstAllocators tstMat tstSP tstAssert tstSortPermutation \
            tstMat1RA tstMat2RA tstMat3RA tstMat4RA tstMat5RA \
	    tstCheck_Strings tstPacking_Utils tstSoft_Equiv, scalar)

AC_DRACO_ENV

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/ds++ configure.in
dnl-------------------------------------------------------------------------dnl
