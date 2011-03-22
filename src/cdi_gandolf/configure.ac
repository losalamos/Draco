dnl-------------------------------------------------------------------------dnl
dnl draco/src/cdi_gandolf configure.ac
dnl Kelly Thompson
dnl Thu Jul 13 16:18:14 2000
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for cdi_gandolf
dnl-------------------------------------------------------------------------dnl

dnl Everthing between AC_INIT and AC_OUTPUT will be used to created
dnl the files listed on the AC_OUTPUT line.

dnl Process any command line arguments and find the source code
dnl directory.
AC_INIT

dnl As an argument to this macro provide the name of a file that
dnl exists in the local package directory.
AC_CONFIG_SRCDIR([GandolfGrayOpacity.hh])

dnl Tell autoconf the directory that contains user defined autoconf
dnl macros needed by this package.  Also tells autoconf to use the
dnl scripts 'install-sh', 'config.sub' and 'config.guess' found in
dnl this directory.
AC_CONFIG_AUX_DIR(../../config)

# Toggle between autoconf and cmake build systems when constructing
# config.h from config.h.in
AC_DEFINE(AUTOCONF autoconf)

dnl If needed, AC_CONFIG_HEADER(config.h.in) should go here (p. 15)
dnl This will replace @DEFS@ with -DHAVE_CONFIG_H in generated files
dnl and will create the file config.h with a list of #define val1
dnl statements based on AC_DEFINE(CPP_VAR, value) commands.
AC_CONFIG_HEADERS(cdi_gandolf/config.h:config.h.in)
dnl The text string that will identify this package to the build
dnl system. 
package="cdi_gandolf"

dnl VENDOR SETUPS
dnl [KT 2010-4-28] Note that the setup done in AC_ALL_VENDORS_SETUP
dnl    does not propigate to the sub-packages.  Here we register the
dnl    gandolf vendor as required for building the unit tests.
AC_GANDOLF_SETUP(test)dnl vendor is only used for the unit tests.

dnl Append this directory to a list of header install locations.
AC_INSTALL_HEADERS

dnl Append this directory to a list of library install locations.
AC_INSTALL_LIB

dnl A space delimited list of packages that this package depends on.
AC_NEEDS_LIBS(cdi ds++)

dnl DRACO SETUP

dnl scalar tests
AC_RUNTESTS(tGandolfFile tGandolfOpacity tGandolfWithCDI tGandolfOdfmgOpacity, scalar)
AC_RUNTESTS(ReadOdfGandolfFile, app)

dnl DRACO SETUP

dnl Evaluate requested build environment and create values to be
dnl substituted into the build system files listed in AC_OUTPUT
AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-6_1_0)

dnl Create build system files in the ${prefix} directory.
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
