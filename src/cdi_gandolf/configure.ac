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
dnl [KT 2010-4-28] Note that the setup done in AC_ALL_VENDORS_SETUP
dnl    does not propigate to the sub-packages.  Here we register the
dnl    gandolf vendor as required for building the unit tests.
AC_GANDOLF_SETUP(test)dnl vendor is only used for the unit tests.

AC_INSTALL_HEADERS
AC_INSTALL_LIB
AC_NEEDS_LIBS(cdi ds++)

AC_RUNTESTS(tGandolfFile tGandolfOpacity tGandolfWithCDI tGandolfOdfmgOpacity, scalar)
AC_RUNTESTS(ReadOdfGandolfFile, app)

dnl DRACO SETUP
AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-5_20_0)

AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in \
	         doc/Makefile:../../config/Makefile.doc.in])
AC_OUTPUT

# The IPCRESS files are binary so we copy them to the target directory manually.
echo copying ${srcdir}/test/Al_BeCu.ipcress to test/Al_BeCu.ipcress
cp --backup=simple  ${srcdir}/test/Al_BeCu.ipcress test/Al_BeCu.ipcress

echo copying ${srcdir}/test/analyticOpacities.ipcress to test/analyticOpacities.ipcress
cp --backup=simple ${srcdir}/test/analyticOpacities.ipcress test/analyticOpacities.ipcress

echo copying ${srcdir}/test/odfregression10.ipcress to test/odfregression10.ipcress
cp --backup=simple  ${srcdir}/test/odfregression10.ipcress test/odfregression10.ipcress

# If build directory == srcdir, (in-tree build) cp will make
# backups. Remove them if they exist.
rm -f test/Al_BeCu.ipcress~ test/analyticOpacities.ipcress~ test/odfregression10.ipcress~

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/cdi_gandolf configure.ac
dnl-------------------------------------------------------------------------dnl
