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

AC_RUNTESTS(tstCheck_Strings \
            tstAssert \
            tstShared_Lib \
            tstSoft_Equiv \
            tstSP \
            tstRCF \
            tstPacking_Utils \
            tstAllocators \
            tstSortPermutation \
            tstVector_Lite \
            tstMat \
            tstMat1RA \
            tstMat2RA \
            tstMat3RA \ 
            tstMat4RA \
            tstMat5RA \
            , scalar)

AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-4_3_0)

dnl ---------------------------------------------------
dnl Need libdl for dlopen, which is used by Shared_Lib.
dnl ---------------------------------------------------
#
# Case through different hosts, if it is a host we haven't configured
# yet then make the code "NULL"
#
case $host in

   # >>> LINUX
   *-linux-gnu)
       LIBS="$LIBS -ldl"
   ;;
   
   # >>> TRU64 COMPAQ
   #alpha*-dec-osf*)
       
   #;;

   # >>> All others
   *)
       AC_DEFINE(NO_DLOPEN)
   ;;
esac

dnl -------------------------------------------
dnl look for flexelint, dot and ps2pdf on Linux
dnl -------------------------------------------
AC_PATH_PROG(FLEXELINT_PATH, lin, null, /usr/local/flexelint/bin )
AC_PATH_PROG(DOTCMD, dot, null)
AC_PATH_PROG(PS2PDF, ps2pdf, null)

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in \
                 doc/Makefile:../../config/Makefile.doc.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/ds++ configure.in
dnl-------------------------------------------------------------------------dnl
