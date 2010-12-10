dnl-------------------------------------------------------------------------dnl
dnl File draco/src/special_functions configure.ac
dnl Kelly Thompson
dnl Tue Feb 22 09:48:39 2000
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.ac for special_functions
dnl-------------------------------------------------------------------------dnl

dnl Everthing between AC_INIT and AC_OUTPUT will be used to created
dnl the files listed on the AC_OUTPUT line.

dnl Process any command line arguments and find the source code
dnl directory.
AC_INIT
AC_CONFIG_SRCDIR([Release.hh])

dnl echo " "
dnl echo "======================================================================"
dnl echo "Running configure for draco/src/special_functions"
dnl echo "======================================================================"

dnl Tell autoconf the directory that contains user defined autoconf
dnl macros needed by this package.  Also tells autoconf to use the
dnl scripts 'install-sh', 'config.sub' and 'config.guess' found in
dnl this directory.
AC_CONFIG_AUX_DIR(../../config)

dnl If needed, AC_CONFIG_HEADER(config.h.in) should go here (p. 15)
dnl This will replace @DEFS@ with -DHAVE_CONFIG_H in generated files
dnl and will create the file config.h with a list of #define val1
dnl statements based on AC_DEFINE(CPP_VAR, value) commands.
AC_CONFIG_HEADER(special_functions/config.h:config.h.in)

dnl The text string that will identify this package to the build
dnl system. 
package="special_functions"

dnl
dnl VENDOR DEPENDENCIES
dnl 
dnl Register vendors dependencies for this component
AC_LAPACK_SETUP(pkg)
AC_GSL_SETUP(pkg)

dnl Append this directory to a list of header install locations.
AC_INSTALL_HEADERS

dnl Append this directory to a list of library install locations.
AC_INSTALL_LIB

dnl A space delimited list of packages that this package depends on.
AC_NEEDS_LIBS( ds++ units ode )

dnl tests
AC_RUNTESTS(\
   tstSpecialFunctions \
   tstYlm \
   tstPower \
   tstFM12 \
   tstF12inv \
   tstF12 \
   tstF1 \
   tstF2inv \
   tstF32 \
   tstF3 \
   tstF4 \
   tstF_eta_inv \
   , scalar)

dnl Evaluate requested build environment and create values to be
dnl substituted into the build system files listed in AC_OUTPUT
AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-6_0_0)

dnl Create build system files in the ${prefix} directory.
AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl		 end of draco/src/special_functions configure.ac
dnl-------------------------------------------------------------------------dnl
