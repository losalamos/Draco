dnl-------------------------------------------------------------------------dnl
dnl draco/src/cdi_gandolf configure.ac
dnl Kelly Thompson
dnl Thu Jul 13 16:18:14 2000
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for cdi_gandolf
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([GandolfGrayOpacity.hh])
AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_HEADERS(cdi_gandolf/config.h:config.h.in)
package="cdi_gandolf"

dnl VENDOR SETUPS
dnl why again - already done via draco/configure.in AC_ALL_VENDORS_SETUP
AC_GANDOLF_SETUP(pkg)

AC_INSTALL_HEADERS
AC_INSTALL_LIB
AC_NEEDS_LIBS(cdi ds++)

AC_RUNTESTS(tGandolfFile tGandolfOpacity tGandolfWithCDI, scalar)

dnl DRACO SETUP
AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-5_16_0)

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in \
	         doc/Makefile:../../config/Makefile.doc.in])
AC_OUTPUT

# The IPCRESS files are binary so we copy them to the target directory manually.
echo copying ${srcdir}/test/Al_BeCu.ipcress to test/Al_BeCu.ipcress
cp -f --backup=simple  ${srcdir}/test/Al_BeCu.ipcress test/Al_BeCu.ipcress
echo copying ${srcdir}/test/analyticOpacities.ipcress to test/analyticOpacities.ipcress
cp -f --backup=simple ${srcdir}/test/analyticOpacities.ipcress test/analyticOpacities.ipcress

# If build directory == srcdir, (in-tree build) cp will make
# backups. Remove them if they exist.
rm -f test/Al_BeCu.ipcres~ test/analyticOpacities.ipcress~

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/cdi_gandolf configure.ac
dnl-------------------------------------------------------------------------dnl
