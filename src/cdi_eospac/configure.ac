dnl-------------------------------------------------------------------------dnl
dnl draco/src/cdi_eospac configure.in
dnl Kelly Thompson
dnl Fri Mar 30 10:43:31 2001
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for cdi_eospac
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([Eospac.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER([cdi_eospac/config.h:config.h.in])
package="cdi_eospac"

dnl VENDOR SETUPS
dnl Setup up Eospac.  The package needs the libs.
AC_EOSPAC5_SETUP(pkg)

AC_INSTALL_HEADERS
AC_INSTALL_LIB
AC_NEEDS_LIBS(cdi ds++)

AC_RUNTESTS(tEospac tEospacWithCDI, scalar)

dnl DRACO SETUP
AC_DRACO_ENV

dnl Setup config.h for EOSPAC type defines.
case $host in
mips-sgi-irix6.*)
   AC_DEFINE(V_INT, int)
   AC_DEFINE(V_FLOAT, double)
   AC_DEFINE(V_BOOL, int)
   ;;
*-pc-linux-gnu)
   AC_DEFINE(V_INT, long)
   AC_DEFINE(V_FLOAT, double)
   AC_DEFINE(V_BOOL, long) 
   ;;
alpha*-dec-osf*)
   AC_DEFINE(V_INT, int)
   AC_DEFINE(V_FLOAT, double)
   AC_DEFINE(V_BOOL, int)
esac

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in \
                 doc/Makefile:../../config/Makefile.doc.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl end of draco/src/cdi_eospac configure.in
dnl-------------------------------------------------------------------------dnl
