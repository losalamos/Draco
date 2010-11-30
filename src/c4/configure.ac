dnl-------------------------------------------------------------------------dnl
dnl c4/configure.ac
dnl Thomas M. Evans
dnl Fri Jan  8 14:46:27 1999
dnl-------------------------------------------------------------------------dnl
dnl $Id$
dnl-------------------------------------------------------------------------dnl

dnl Everthing between AC_INIT and AC_OUTPUT will be used to created
dnl the files listed on the AC_OUTPUT line.

dnl Process any command line arguments and find the source code
dnl directory.
AC_INIT

dnl As an argument to this macro provide the name of a file that
dnl exists in the local package directory. 
AC_CONFIG_SRCDIR([C4_Req.hh])

dnl echo " "
dnl echo "======================================================================"
dnl echo "Running configure for draco/src/c4"
dnl echo "======================================================================"

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
AC_CONFIG_HEADER(c4/config.h:config.h.in)
dnl The text string that will identify this package to the build
dnl system. 
package="c4"

dnl
dnl VENDOR DEPENDENCIES
dnl 
dnl Register vendors dependencies for this component
AC_MPI_SETUP(pkg)

dnl Append this directory to a list of header install locations.
AC_INSTALL_HEADERS

dnl Append this directory to a list of library install locations.
AC_INSTALL_LIB

dnl A space delimited list of packages that this package depends on.
AC_NEEDS_LIBS(ds++)

dnl DRACO SETUP

dnl parallel tests
AC_RUNTESTS( \
            tstC4_Req \
            tstTime \
            tstReduction \
            tstPingPong \
            tstComm_Dup \
            tstBroadcast \
            tstAbort \
            tstScalar \
            tstParallelUnitTest \
            tstInvert_Comm_Map \
            tstCompare \
            tstSend_Receive \
            tstSwap \
            tstGatherScatter \
            tstProcessor_Group \
            tstTermination_Detector \
            tstglobal_containers \
            , 1 2 4)

dnl application tests
AC_TEST_APPLICATION( tstApplicationUnitTest, 1 2 4 )

dnl Determine system type.  These defs are saved in utils/config.h.
case `uname` in
Linux)
  AC_DEFINE(c4_isLinux) ;; 
OSF1)
  AC_DEFINE(c4_isOSF1) ;; 
Darwin)
  AC_DEFINE(c4_isDarwin) ;;
AIX)
  AC_DEFINE(c4_isAIX) ;;
esac

dnl Evaluate requested build environment and create values to be
dnl substituted into the build system files listed in AC_OUTPUT
AC_DRACO_ENV
AC_DRACO_AUTODOC(draco-5_21_0)

dnl Create build system files in the ${prefix} directory.
AC_CONFIG_FILES([Makefile:../../config/Makefile.package.in \
                 test/Makefile:../../config/Makefile.test.in])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of c4/configure.ac
dnl-------------------------------------------------------------------------dnl
