dnl-------------------------------------------------------------------------dnl
dnl draco/src/cdi_eospac configure.ac
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
AC_DRACO_AUTODOC

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

   dnl Lahey 6.2 and 6.0 are not compatible.  EOSPAC 5.35 libraries
   dnl are built with 6.0.  If we are not using Lahey 6.2 provide a 
   dnl configure time error.

dnl   AC_MSG_CHECKING("for EOSPAC v5.35 compatible version of F90")
dnl   cdi_eospac_lf95=`which lf95`
dnl   if test -z `echo ${cdi_eospac_lf95} | grep lf9562`; then
dnl      AC_MSG_RESULT("ok. found ${cdi_eospac_lf95}")
dnl   else
dnl      AC_MSG_RESULT("Lahey 6.2 (${cdi_eospac_lf95}) is not compatible with EOSPAC v5.35")
dnl      AC_MSG_ERROR("EOSPAC v5.35 is incompatible with ${cdi_eospac_lf95}.  Use Lahey 6.0 or 6.1 instead.")
dnl   fi

   dnl Even if we are using Lahey 6.2, cdi_eospac will work if we link
   dnl against the Lahey 6.0 libraries.
   dnl RPATH="${RPATH} -Xlinker -rpath /usr/local/lf9560/lib"
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
