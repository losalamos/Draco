dnl-------------------------------------------------------------------------dnl
dnl draco/src/plot2D configure.in
dnl Rob Lowrie
dnl Mon Apr 15 10:05:15 2002
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for plot2D
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([Plot2D.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADER(plot2D/config.h:config.h.in)
package="plot2D"

dnl VENDOR SETUPS
AC_GRACE_SETUP(pkg)

AC_INSTALL_HEADERS
AC_INSTALL_LIB
AC_NEEDS_LIBS(ds++)
 
dnl tests
AC_RUNTESTS(tstPlot2D, scalar)

dnl DRACO SETUP
AC_DRACO_ENV

dnl Grace must be in the path!
AC_PATH_PROG(XMGRACE_PATH, xmgrace, null)
if test "${XMGRACE_PATH}" = null; then
   AC_MSG_ERROR("No valid xmgrace found!")
fi

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in \
                 test/tstPlot2D_Diff.py:test/python/tstPlot2D_Diff.py.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/plot2D configure.in
dnl-------------------------------------------------------------------------dnl
