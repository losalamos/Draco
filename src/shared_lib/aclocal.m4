# generated automatically by aclocal 1.7.3 -*- Autoconf -*-

# Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002
# Free Software Foundation, Inc.
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

dnl-------------------------------------------------------------------------dnl
dnl ac_conf.m4
dnl
dnl Service macros used in configure.ac's throughout Draco.
dnl
dnl Thomas M. Evans
dnl 1999/02/04 01:56:19
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl AC_DRACO_PREREQ
dnl
dnl Checks the configure version
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DRACO_PREREQ], [dnl

   # we need at least autoconf 2.53 to work correctly
   AC_PREREQ(2.53)

])

dnl-------------------------------------------------------------------------dnl
dnl AC_NEEDS_LIBS
dnl
dnl add DRACO-dependent libraries necessary for a package
dnl usage: configure.ac
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_NEEDS_LIBS, [dnl
   if test ${has_libdir:=no} != "yes" ; then
       DRACO_LIBS="${DRACO_LIBS} -L\${libdir}"
       has_libdir="yes"
   fi

   for lib in $1
   do
       # temporary string to keep line from getting too long
       draco_depends="\${libdir}/lib\${LIB_PREFIX}${lib}\${libsuffix}"
       DRACO_DEPENDS="${DRACO_DEPENDS} ${draco_depends}"
       DRACO_LIBS="${DRACO_LIBS} -l\${LIB_PREFIX}${lib}"
   done

   # Keep a list of component dependencies free of other tags or paths.
   DEPENDENT_COMPONENTS="$1"

])

dnl-------------------------------------------------------------------------dnl
dnl AC_NEEDS_LIBS_TEST
dnl
dnl add DRACO-dependent libraries necessary for a package test
dnl usage: configure.ac
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_NEEDS_LIBS_TEST, [dnl
   DRACO_TEST_LIBS="${DRACO_TEST_LIBS} -L\${libdir}"
   for lib in $1
   do
       # temporary string to keep line from getting too long
       draco_test_depends="\${libdir}/lib\${LIB_PREFIX}${lib}\${libsuffix}"
       DRACO_TEST_DEPENDS="${DRACO_TEST_DEPENDS} ${draco_test_depends}"
       DRACO_TEST_LIBS="${DRACO_TEST_LIBS} -l\${LIB_PREFIX}${lib}"
   done
])

dnl-------------------------------------------------------------------------dnl
dnl AC_RUNTESTS
dnl
dnl add DRACO-package tests (default to use DejaGnu)
dnl usage: in configure.ac:
dnl AC_RUNTESTS(testexec1 testexec2 ... , {nprocs1 nprocs2 ... | scalar})
dnl where serial means run as serial test only.
dnl If compiling with scalar c4 then nprocs are ignored.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_RUNTESTS, [dnl
	test_alltarget="$test_alltarget $1"
        
	test_nprocs="$2"

	if test -z "${test_nprocs}" ; then
	    AC_MSG_ERROR("No procs choosen for the tests!")
        fi
])

dnl-------------------------------------------------------------------------dnl
dnl AC_TESTEXE
dnl
dnl determines what type of executable the tests are, for example, you 
dnl can set the executable to some scripting extension, like python.
dnl the default is an executable binary
dnl options are PYTHON
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_TESTEXE, [dnl
   test_exe="$1"
])

dnl-------------------------------------------------------------------------dnl
dnl AC_INSTALL_EXECUTABLE
dnl
dnl where executables will be installed
dnl usage: configure.ac
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_INSTALL_EXECUTABLE, [ dnl
   install_executable="\${bindir}/\${package}"
   installdirs="${installdirs} \${bindir}"
   alltarget="${alltarget} bin/\${package}"
])

dnl-------------------------------------------------------------------------dnl
dnl AC_INSTALL_LIB
dnl
dnl where libraries will be installed
dnl usage: configure.ac
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_INSTALL_LIB, [ dnl
   install_lib='${libdir}/lib${LIB_PREFIX}${package}${libsuffix}'
   installdirs="${installdirs} \${libdir}"
   alltarget="${alltarget} lib\${LIB_PREFIX}\${package}\${libsuffix}"

   # test will need to link this library
   PKG_DEPENDS='../lib${LIB_PREFIX}${package}${libsuffix}'
   PKG_LIBS='-L.. -l${LIB_PREFIX}${package}'
])

dnl-------------------------------------------------------------------------dnl
dnl AC_INSTALL_HEADERS
dnl
dnl where headers will be installed 
dnl usage: configure.ac
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_INSTALL_HEADERS, [ dnl
   install_headers="\${installheaders}"
   installdirs="${installdirs} \${includedir} \${includedir}/\${package}"
])

dnl-------------------------------------------------------------------------dnl
dnl AC_CHECK_TOOLS
dnl
dnl Find tools used by the build system (latex, bibtex, python, etc)
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DRACO_CHECK_TOOLS], [dnl

   dnl
   dnl TOOL CHECKS
   dnl
   
   dnl check for and assign the path to python
   AC_PATH_PROG(PYTHON_PATH, python, null)
   if test "${PYTHON_PATH}" = null ; then
       AC_MSG_ERROR("No valid Python found!")
   fi
   
   dnl check for and assign the path to perl
   AC_PATH_PROG(PERL_PATH, perl, null)
   if test "${PERL_PATH}" = null ; then
       AC_MSG_WARN("No valid Perl found!")
   fi

   dnl check for CVS
   AC_PATH_PROG(CVS_PATH, cvs, null)
   if test "${CVS_PATH}" = null ; then
       AC_MSG_WARN("No valid CVS found!")
   fi

   dnl check for and assign the path to ghostview
   AC_CHECK_PROGS(GHOSTVIEW, ghostview gv, null)
   if test "${GHOSTVIEW}" = null ; then
       AC_MSG_WARN("No valid ghostview found!")
   fi

   dnl check for and assign the path to latex
   AC_CHECK_PROGS(LATEX, latex, null)
   if test "${LATEX}" = null ; then
       AC_MSG_WARN("No valid latex found!")
   fi
   AC_SUBST(LATEXFLAGS)

   dnl check for and assign the path to bibtex
   AC_CHECK_PROGS(BIBTEX, bibtex, null)
   if test "${BIBTEX}" = null ; then
       AC_MSG_WARN("No valid bibtex found!")
   fi
   AC_SUBST(BIBTEXFLAGS)

   dnl check for and assign the path to xdvi
   AC_CHECK_PROGS(XDVI, xdvi, null)
   if test "${XDVI}" = null ; then
       AC_MSG_WARN("No valid xdvi found!")
   fi
   AC_SUBST(XDVIFLAGS)

   dnl check for and assign the path to dvips
   AC_CHECK_PROGS(DVIPS, dvips, null)
   if test "${DVIPS}" = null ; then
       AC_MSG_WARN("No valid dvips found!")
   fi
   AC_SUBST(DVIPSFLAGS)

   dnl check for and assign the path for printing (lp)
   AC_CHECK_PROGS(LP, lp lpr, null)
   if test "${LP}" = null ; then
       AC_MSG_WARN("No valid lp or lpr found!")
   fi
   AC_SUBST(LPFLAGS)

   dnl check for and assign the path for doxygen
   AC_PATH_PROG(DOXYGEN_PATH, doxygen, null)
   if test "${DOXYGEN_PATH}" = null ; then
       AC_MSG_WARN("No valid Doxygen found!")
   fi
   AC_SUBST(DOXYGEN_PATH)

])

dnl-------------------------------------------------------------------------dnl
dnl AC_ASCI_WHITE_TEST_WORK_AROUND_PREPEND
dnl
dnl changes compiler from newmpxlC to newxlC so that tests can be run
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_ASCI_WHITE_TEST_WORK_AROUND_PREPEND], [dnl

   # change compiler
   if test "${CXX}" = newmpxlC; then
       white_compiler='newmpxlC'
       CXX='newxlC'
       AC_MSG_WARN("Changing to ${CXX} compiler for configure tests.")
   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_ASCI_WHITE_TEST_WORK_AROUND_APPEND
dnl
dnl changes compiler back to newmpxlC
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_ASCI_WHITE_TEST_WORK_AROUND_APPEND], [dnl

   # change compiler back
   if test "${white_compiler}" = newmpxlC; then
       CXX='newmpxlC'
       AC_MSG_WARN("Changing back to ${CXX} compiler.")
   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_HEAD_MAKEFILE
dnl 
dnl Builds default makefile in the head directory
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_HEAD_MAKEFILE], [dnl

   AC_FIND_TOP_SRC($srcdir, package_top_srcdir)
   AC_DBS_VAR_SUBSTITUTIONS
   AC_CONFIG_FILES([Makefile:config/Makefile.head.in])

])

dnl-------------------------------------------------------------------------dnl
dnl AC_SRC_MAKEFILE
dnl 
dnl Builds default makefile in the src directory
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_SRC_MAKEFILE], [dnl

   AC_FIND_TOP_SRC($srcdir, package_top_srcdir)
   AC_DBS_VAR_SUBSTITUTIONS
   AC_CONFIG_FILES([Makefile:../config/Makefile.src.in])

])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_conf.m4
dnl-------------------------------------------------------------------------dnl


dnl-------------------------------------------------------------------------dnl
dnl ac_local.m4
dnl
dnl Macros used internally within the Draco build system.
dnl
dnl Thomas M. Evans
dnl 1999/02/04 01:56:22
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl AC_WITH_DIR
dnl
dnl Define --with-xxx[=DIR] with defaults to an environment variable.
dnl       Usage: AC_WITH_DIR(flag, CPPtoken, DefaultValue, HelpStr)
dnl                for environment variables enter \${ENVIRONVAR} for
dnl                DefaultValue
dnl usage: in aclocal.m4
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_WITH_DIR, [dnl

 dnl
 dnl  The following M4 macros will be expanded into the body of AC_ARG_WITH
 dnl
 dnl AC_PACKAGE is the flag with all dashes turned to underscores
 dnl AC_WITH_PACKAGE will be substituted to the autoconf shell variable
 dnl    with_xxx
 dnl AC_CMDLINE is the shell command to strip double and trailing slashes
 dnl    from directory names.

 define([AC_PACKAGE], [translit($1, [-], [_])])dnl
 define([AC_WITH_PACKAGE], [with_]AC_PACKAGE)dnl
 define([AC_CMDLINE],dnl
[echo "$]AC_WITH_PACKAGE[" | sed 's%//*%/%g' | sed 's%/$%%'])dnl

 AC_ARG_WITH($1,
   [  --with-$1[=DIR]    $4 ($3 by default)],
   if test $AC_WITH_PACKAGE != "no" ; then
      if test $AC_WITH_PACKAGE = "yes" ; then
         # following eval needed to remove possible '\' from $3
         eval AC_WITH_PACKAGE=$3
      fi

      # this command removes double slashes and any trailing slash

      AC_WITH_PACKAGE=`eval AC_CMDLINE`
      if test "$AC_WITH_PACKAGE:-null}" = "null" ; then
         { echo "configure: error: --with-$1 directory is unset" 1>&2; \
           exit 1; }
      fi
      if test ! -d $AC_WITH_PACKAGE ; then
         { echo "configure: error: $AC_WITH_PACKAGE: invalid directory" 1>&2; \
           exit 1; }
      fi

      # this sets up the shell variable, with the name of the CPPtoken,
      # and that we later will do an AC_SUBST on.
      $2="${AC_WITH_PACKAGE}/"

      # this defines the CPP macro with the directory and single slash appended.
      AC_DEFINE_UNQUOTED($2, ${AC_WITH_PACKAGE}/)dnl

      # print a message to the users (that can be turned off with --silent)

      echo "$2 has been set to $$2" 1>&6

   fi)

   AC_SUBST($2)dnl

])
	
dnl-------------------------------------------------------------------------dnl
dnl AC_VENDORLIB_SETUP(1,2)
dnl
dnl set up for VENDOR_LIBS or VENDOR_TEST_LIBS
dnl usage: in aclocal.m4
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_VENDORLIB_SETUP, [dnl

   # $1 is the vendor_<> tag (equals pkg or test)
   # $2 are the directories added 

   if test "${$1}" = pkg ; then
       VENDOR_LIBS="${VENDOR_LIBS} $2"
   elif test "${$1}" = test ; then
       VENDOR_TEST_LIBS="${VENDOR_TEST_LIBS} $2"
   fi
])

dnl-------------------------------------------------------------------------dnl
dnl AC_FIND_TOP_SRC(1,2)
dnl 
dnl Find the top source directory of the package by searching upward
dnl from the argument directory. The top source directory is defined
dnl as the one with a 'config' sub-directory.
dnl
dnl Note: This function will eventually quit if the searched for
dnl directory is not above the argument. It does so when $temp_dir
dnl ceases to be a valid directory, which only seems to happen after a
dnl LOT of ..'s are added to it.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_FIND_TOP_SRC, [dnl
   
   # $1 is the component's source directory
   # $2 is the variable to store the package's main source directory in.

   temp_dir=$1
   AC_MSG_CHECKING([package top source directory])
   while test -d $temp_dir -a ! -d $temp_dir/config ; do   
       temp_dir="${temp_dir}/.."
   done
   if test -d $temp_dir; then
       $2=`cd $temp_dir; pwd;`
       AC_MSG_RESULT([$$2])
   else
       AC_MSG_ERROR('Could not find package top source directory')
   fi
])

dnl-------------------------------------------------------------------------dnl
dnl DO VARIABLE SUBSTITUTIONS ON AC_OUTPUT
dnl
dnl These are all the variable substitutions used within the draco
dnl build system
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_VAR_SUBSTITUTIONS], [dnl

   # these variables are declared "precious", meaning that they are
   # automatically substituted, put in the configure --help, and
   # cached 
   AC_ARG_VAR(CC)dnl
   AC_ARG_VAR(CFLAGS)dnl

   AC_ARG_VAR(CXX)dnl
   AC_ARG_VAR(CXXFLAGS)dnl

   AC_ARG_VAR(LD)dnl
   AC_ARG_VAR(LDFLAGS)dnl

   AC_ARG_VAR(AR)dnl
   AC_ARG_VAR(ARFLAGS)dnl

   AC_ARG_VAR(CPPFLAGS)dnl

   # dependency rules
   AC_SUBST(DEPENDENCY_RULES)

   # other compiler substitutions
   AC_SUBST(STRICTFLAG)dnl
   AC_SUBST(PARALLEL_FLAG)dnl
   AC_SUBST(RPATH)dnl
   AC_SUBST(LIB_PREFIX)dnl

   # install program
   AC_SUBST(INSTALL)dnl
   AC_SUBST(INSTALL_DATA)dnl

   # files to install
   : ${installfiles:='${install_executable} ${install_lib} ${install_headers}'}
   AC_SUBST(installfiles)dnl
   AC_SUBST(install_executable)dnl
   AC_SUBST(install_lib)dnl
   AC_SUBST(install_headers)dnl
   AC_SUBST(installdirs)dnl

   # package libraries
   AC_SUBST(alltarget)dnl
   AC_SUBST(libsuffix)dnl
   AC_SUBST(dirstoclean)dnl
   AC_SUBST(package)dnl
   AC_SUBST(DRACO_DEPENDS)dnl
   AC_SUBST(DRACO_LIBS)dnl
   AC_SUBST(VENDOR_DEPENDS)dnl
   AC_SUBST(VENDOR_INC)dnl
   AC_SUBST(VENDOR_LIBS)dnl
   AC_SUBST(ARLIBS)dnl

   # package testing libraries
   AC_SUBST(PKG_DEPENDS)dnl
   AC_SUBST(PKG_LIBS)dnl
   AC_SUBST(DRACO_TEST_DEPENDS)dnl
   AC_SUBST(DRACO_TEST_LIBS)dnl
   AC_SUBST(VENDOR_TEST_DEPENDS)dnl
   AC_SUBST(VENDOR_TEST_LIBS)dnl
   AC_SUBST(ARTESTLIBS)dnl
   AC_SUBST(test_alltarget)dnl
   AC_SUBST(test_flags)dnl
   AC_SUBST(test_scalar)dnl
   AC_SUBST(test_nprocs)dnl
   AC_SUBST(test_output_files)dnl

   # libraries
   AC_ARG_VAR(LIBS)dnl

   # configure options
   AC_SUBST(configure_command)dnl

   # directories in source tree
   AC_SUBST(package_top_srcdir)
   
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_local.m4
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl ac_vendors.m4
dnl
dnl Macros for each vendor that is used supported by the Draco build
dnl system.
dnl
dnl Thomas M. Evans
dnl 1999/02/04 01:56:22
dnl-------------------------------------------------------------------------dnl


dnl-------------------------------------------------------------------------dnl
dnl AC_MPI_SETUP
dnl
dnl MPI implementation (off by default)
dnl MPI is an optional vendor
dnl
dnl we wait to set the basic MPI libraries (if it is on) until
dnl after checking the C4 status; these functions are performed
dnl in ac_dracoenv.m4, section SYSTEM-SPECIFIC SETUP; we do this
dnl here because each platform has different mpi options for
dnl vendors and mpich
dnl
dnl note that we used to do this in a function called AC_COMM_SET;
dnl however, there are too many platform-dependent variables 
dnl to continue doing this; so we do all these operations in the
dnl platform specific section of ac_dracoenv.m4
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_MPI_SETUP], [dnl

   dnl define --with-mpi
   AC_ARG_WITH(mpi,
      [  --with-mpi=[vendor,mpich,lampi] 
	                  determine MPI implementation (vendor on SGI,SUN; mpich on LINUX)])

   dnl define --with-mpi-inc and --with-mpi-lib
   AC_WITH_DIR(mpi-inc, MPI_INC, \${MPI_INC_DIR},
	       [tell where MPI includes are])
   AC_WITH_DIR(mpi-lib, MPI_LIB, \${MPI_LIB_DIR},
	       [tell where MPI libs are])

   # determine if this package is needed for testing or for the
   # package
   vendor_mpi=$1 

   # set default value for with_mpi which is no
   if test "${with_mpi:=no}" = yes ; then 
       with_mpi='vendor'
   fi

   # if the user sets MPI_INC and MPI_LIB directories then turn on  
   # with_mpi and set it to vendor if with_mpi=no to begin with
   if test "${with_mpi}" = no ; then
       if test -n "${MPI_INC}" ; then
	   with_mpi='vendor'
       elif test -n "${MPI_LIB}" ; then
	   with_mpi='vendor'
       fi
   fi
   
   # if c4=mpi and with-mpi=no explicitly then 
   # define them (mpi gets set to vendor by default)
   if test "$with_c4" = mpi ; then
       if test "$with_mpi" = no ; then
	   with_mpi='vendor'
       fi
   fi

]) 


AC_DEFUN([AC_MPI_FINALIZE], [dnl

   # only add stuff if mpi is not no and the vendor is defined
   if test "${with_mpi}" != no && test -n "${vendor_mpi}"; then

       # include path
       if test -n "${MPI_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${MPI_INC}"
       fi
   
       # libraries
       if test -n "${MPI_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_mpi, -L${MPI_LIB} ${mpi_libs})
       elif test -z "${MPI_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_mpi, ${mpi_libs})
       fi

       # add MPI directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${MPI_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${MPI_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_SPRNG_SETUP
dnl
dnl SPRNG LIBRARY SETUP (on by default -lfg)
dnl SPRNG is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_SPRNG_SETUP], [dnl

   dnl define --with-sprng
   AC_ARG_WITH(sprng,
      [  --with-sprng[=lib]      determine the rng lib (lfg is default)])
	
   dnl define --with-sprng-inc and --with-sprng-lib
   AC_WITH_DIR(sprng-inc, SPRNG_INC, \${SPRNG_INC_DIR},
	       [tell where SPRNG includes are])
   AC_WITH_DIR(sprng-lib, SPRNG_LIB, \${SPRNG_LIB_DIR},
	       [tell where SPRNG libraries are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_sprng=$1

   # choices are with_sprng = lfg, lcg, yes, or no

   # default (sprng is set to lfg by default)
   if test "${with_sprng:=lfg}" = yes ; then
       with_sprng='lfg'
   fi

])


AC_DEFUN([AC_SPRNG_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_sprng}"; then

       # include path
       if test -n "${SPRNG_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${SPRNG_INC}"
       fi
   
       # libraries
       if test -n "${SPRNG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_sprng, -L${SPRNG_LIB} -l${with_sprng})
       elif test -z "${SPRNG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_sprng, -l${with_sprng})
       fi

       # add sprng directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${SPRNG_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${SPRNG_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_AZTEC_SETUP
dnl
dnl AZTEC SETUP (on by default)
dnl AZTEC is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_AZTEC_SETUP], [dnl

   dnl define --with-aztec
   AC_ARG_WITH(aztec,
      [  --with-aztec=[lib]      determine the aztec lib (aztec is the default)])
 
   dnl define --with-aztec-inc
   AC_WITH_DIR(aztec-inc, AZTEC_INC, \${AZTEC_INC_DIR},
	       [tell where AZTEC includes are])

   dnl define --with-aztec-lib
   AC_WITH_DIR(aztec-lib, AZTEC_LIB, \${AZTEC_LIB_DIR},
	       [tell where AZTEC libraries are])

   # set default value of aztec includes and libs
   if test "${with_aztec:=aztec}" = yes ; then
       with_aztec='aztec'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_aztec=$1

])


AC_DEFUN([AC_AZTEC_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_aztec}" ; then

       # include path
       if test -n "${AZTEC_INC}"; then 
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${AZTEC_INC}"
       fi

       # library path
       if test -n "${AZTEC_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_aztec, -L${AZTEC_LIB} -l${with_aztec})
       elif test -z "${AZTEC_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_aztec, -l${with_aztec})
       fi

       # add AZTEC directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${AZTEC_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${AZTEC_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_GSL_SETUP
dnl
dnl GSL SETUP (on by default)
dnl GSL is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_GSL_SETUP, [dnl

   dnl define --with-gsl
   AC_ARG_WITH(gsl,
      [  --with-gsl=[gsl] 
                       determine GSL lib (gsl is default)])
 
   dnl define --with-gsl-inc
   AC_WITH_DIR(gsl-inc, GSL_INC, \${GSL_INC_DIR},
	       [tell where GSL includes are])

   dnl define --with-gsl-lib
   AC_WITH_DIR(gsl-lib, GSL_LIB, \${GSL_LIB_DIR},
	       [tell where GSL libraries are])

   # set default value of gsl includes and libs
   if test "${with_gsl:=gsl}" = yes ; then
       with_gsl='gsl'
   fi

   # if atlas is available use it's version of cblas, 
   # otherwise use the version provided by GSL
   if test "${with_lapack}" = atlas; then
       gsl_libs='-lgsl'
   else
       gsl_libs='-lgsl -lgslcblas'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_gsl=$1
])


AC_DEFUN([AC_GSL_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_gsl}"; then

       # include path
       if test -n "${GSL_INC}"; then 
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${GSL_INC}"
       fi

       # library path
       if test -n "${GSL_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gsl, -L${GSL_LIB} ${gsl_libs})
       elif test -z "${GSL_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gsl, ${gsl_libs})
       fi

       # add GSL directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${GSL_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${GSL_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_TRILINOS_SETUP
dnl
dnl TRILINOS SETUP (on by default)
dnl TRILINOS is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_TRILINOS_SETUP], [dnl

   dnl define --with-trilinos
   AC_ARG_WITH(trilinos,
      [  --with-trilinos=[lib]    determine the trilinos implementation (aztecoo is default)])
 
   dnl define --with-trilinos-inc
   AC_WITH_DIR(trilinos-inc, TRILINOS_INC, \${TRILINOS_INC_DIR},
	       [tell where TRILINOS includes are])

   dnl define --with-trilinos-lib
   AC_WITH_DIR(trilinos-lib, TRILINOS_LIB, \${TRILINOS_LIB_DIR},
	       [tell where TRILINOS libraries are])

   # set default value of trilinos includes and libs
   if test "${with_trilinos:=aztecoo}" = yes ; then
       with_trilinos='aztecoo'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_trilinos=$1

])


AC_DEFUN([AC_TRILINOS_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_trilinos}" ; then

       # include path
       if test -n "${TRILINOS_INC}"; then 
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${TRILINOS_INC}"
       fi

       # library path
       if test -n "${TRILINOS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_trilinos, -L${TRILINOS_LIB} -l${with_trilinos} -lepetra -lepetraext -lteuchos -ltriutils)
       elif test -z "${TRILINOS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_trilinos, -l${with_trilinos} -lepetra -lepetraext -lteuchos -ltriutils)
       fi

       # add TRILINOS directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${TRILINOS_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${TRILINOS_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_SCALAPACK_SETUP
dnl
dnl SCALAPACK SETUP 
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_SCALAPACK_SETUP], [dnl

   dnl define --with-scalapack
   AC_ARG_WITH(scalapack,
      [  --with-scalapack=[scalapack] ])
 
   dnl define --with-scalapack-lib
   AC_WITH_DIR(scalapack-lib, SCALAPACK_LIB, \${SCALAPACK_LIB_DIR},
	       [tell where SCALAPACK libraries are])

   # set default value of scalapack includes and libs
   if test "${with_scalapack:=scalapack}" = yes ; then
       with_scalapack='scalapack'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_scalapack=$1

])


AC_DEFUN([AC_SCALAPACK_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_scalapack}" ; then

       # no includes for scalapack

       # library path
       if test -n "${SCALAPACK_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_scalapack, -L${SCALAPACK_LIB} -lscalapack)
       elif test -z "${SCALAPACK_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_scalapack, -lscalapack)
       fi

       # add SCALAPACK directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${SCALAPACK_LIB}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_BLACS_SETUP
dnl
dnl BLACS SETUP 
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_BLACS_SETUP], [dnl

   dnl define --with-blacs
   AC_ARG_WITH(blacs,
      [  --with-blacs=[blacs] ])
 
   dnl define --with-blacs-lib
   AC_WITH_DIR(blacs-lib, BLACS_LIB, \${BLACS_LIB_DIR},
	       [tell where BLACS libraries are])

   # set default value of blacs includes and libs
   if test "${with_blacs:=blacs}" = yes ; then
       with_blacs='blacs'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_blacs=$1

])


AC_DEFUN([AC_BLACS_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_blacs}" ; then

       # no includes for blacs

       # library path
       if test -n "${BLACS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_blacs, -L${BLACS_LIB} -lblacsF77init -lblacsCinit -lblacs -lblacsCinit -lblacs)
       elif test -z "${BLACS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_blacs, -lblacsF77init -lblacsCinit -lblacs -lblacsCinit -lblacs)
       fi

       # add BLACS directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${BLACS_LIB}"

   fi

])
dnl-------------------------------------------------------------------------dnl
dnl AC_HYPRE_SETUP
dnl
dnl HYPRE SETUP 
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_HYPRE_SETUP], [dnl

   dnl define --with-hypre
   AC_ARG_WITH(hypre,
      [  --with-hypre=[hypre] ])
 
   dnl define --with-hypre-inc
   AC_WITH_DIR(hypre-inc, HYPRE_INC, \${HYPRE_INC_DIR},
	       [tell where HYPRE includes are])

   dnl define --with-hypre-lib
   AC_WITH_DIR(hypre-lib, HYPRE_LIB, \${HYPRE_LIB_DIR},
	       [tell where HYPRE libraries are])

   # set default value of hypre includes and libs
   if test "${with_hypre:=hypre}" = yes ; then
       with_hypre='hypre'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_hypre=$1

])


AC_DEFUN([AC_HYPRE_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_hypre}" ; then

       # include path
       if test -n "${HYPRE_INC}"; then 
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${HYPRE_INC}"
       fi

       # library path
       if test -n "${HYPRE_LIB}" ; then

	   AC_VENDORLIB_SETUP(vendor_hypre, -L${HYPRE_LIB} -lHYPRE_parcsr_ls -lHYPRE_DistributedMatrixPilutSolver -lHYPRE_ParaSails -lHYPRE_Euclid -lHYPRE_MatrixMatrix -lHYPRE_DistributedMatrix -lHYPRE_IJ_mv -lHYPRE_parcsr_mv -lHYPRE_seq_mv -lHYPRE_krylov -lHYPRE_utilities)

       elif test -z "${HYPRE_LIB}" ; then

	   AC_VENDORLIB_SETUP(vendor_hypre, -lHYPRE_parcsr_ls -lHYPRE_DistributedMatrixPilutSolver -lHYPRE_ParaSails -lHYPRE_Euclid -lHYPRE_MatrixMatrix -lHYPRE_DistributedMatrix -lHYPRE_IJ_mv -lHYPRE_parcsr_mv -lHYPRE_seq_mv -lHYPRE_krylov -lHYPRE_utilities)

       fi

       # add HYPRE directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${HYPRE_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${HYPRE_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_METIS_SETUP
dnl
dnl METIS SETUP (on by default)
dnl METIS is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_METIS_SETUP], [dnl

   dnl define --with-metis
   AC_ARG_WITH(metis,
      [  --with-metis=[lib]    the metis implementation])
 
   dnl define --with-metis-inc
   AC_WITH_DIR(metis-inc, METIS_INC, \${METIS_INC_DIR},
	       [tell where METIS includes are])

   dnl define --with-metis-lib
   AC_WITH_DIR(metis-lib, METIS_LIB, \${METIS_LIB_DIR},
	       [tell where METIS libraries are])

   # set default value of metis includes and libs
   if test "${with_metis:=metis}" = yes ; then
       with_metis='metis'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_metis=$1

])


AC_DEFUN([AC_METIS_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_metis}" ; then

       # include path
       if test -n "${METIS_INC}"; then 
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${METIS_INC}"
       fi

       # library path
       if test -n "${METIS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_metis, -L${METIS_LIB} -l${with_metis})
       elif test -z "${METIS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_metis, -l${with_metis})
       fi

       # add METIS directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${METIS_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${METIS_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_PCG_SETUP
dnl
dnl PCG LIBRARY SETUP (on by default)
dnl PCG is a required vendor
dnl
dnl note that we add some system-specific libraries for this
dnl vendor in AC_DRACO_ENV; also, the user must set up LAPACK for
dnl this to work
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_PCG_SETUP], [dnl

   dnl define --with-pcg
   AC_ARG_WITH(pcg,        
      [  --with-pcg[=lib]        determine the pcg lib name (pcg is default)])

   dnl define --with-pcg-lib
   AC_WITH_DIR(pcg-lib, PCG_LIB, \${PCG_LIB_DIR},
	       [tell where PCG libraries are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_pcg=$1

   # pcg is set to libpcg by default
   if test "${with_pcg:=pcg}" = yes ; then
       with_pcg='pcg'
   fi

])


AC_DEFUN([AC_PCG_FINALIZE], [dnl

   # set up the libraries
   if test -n "${vendor_pcg}"; then

       # library path
       if test -z "${PCG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_pcg, -l${with_pcg})
       elif test -n "${PCG_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_pcg, -L${PCG_LIB} -l${with_pcg})
       fi

       # add PCG directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${PCG_LIB}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_GANDOLF_SETUP
dnl
dnl GANDOLF LIBRARY SETUP (on by default)
dnl GANDOLF is a required vendor
dnl
dnl SGI needs "-lfortran" on the load line when including libgandolf.a.
dnl This library is added to ${LIBS} in AC_DRACO_ENV.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_GANDOLF_SETUP], [dnl

   dnl define --with-gandolf
   AC_ARG_WITH(gandolf,        
      [  --with-gandolf[=lib]    determine the gandolf lib name (gandolf is default)])

   dnl define --with-gandolf-lib
   AC_WITH_DIR(gandolf-lib, GANDOLF_LIB, \${GANDOLF_LIB_DIR},
	       [tell where GANDOLF libraries are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_gandolf=$1

   # gandolf is set to libgandolf by default
   if test "${with_gandolf:=gandolf}" = yes ; then
       with_gandolf='gandolf'
   fi

])


AC_DEFUN([AC_GANDOLF_FINALIZE], [dnl

   # set up the libraries
   if test -n "${vendor_gandolf}"; then

       # set up library paths
       if test -z "${GANDOLF_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gandolf, -l${with_gandolf})
       elif test -n "${GANDOLF_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gandolf, -L${GANDOLF_LIB} -l${with_gandolf})
       fi

       # add GANDOLF directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${GANDOLF_LIB}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_EOSPAC5_SETUP
dnl
dnl EOSPAC5 LIBRARY SETUP (on by default)
dnl EOSPAC5 is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_EOSPAC5_SETUP], [dnl

   dnl define --with-eospac
   AC_ARG_WITH(eospac,        
      [  --with-eospac[=lib]     determine the eospac lib name (eospac is default)])

   dnl define --with-eospac-lib
   AC_WITH_DIR(eospac-lib, EOSPAC5_LIB, \${EOSPAC5_LIB_DIR},
	       [tell where EOSPAC5 libraries are])

   # determine if this package is needed for testing or for the 
   # package (valid values are pkg or test)
   vendor_eospac=$1

   # eospac is set to libeospac by default
   if test "${with_eospac:=eospac}" = yes ; then
       with_eospac='eospac'
   fi

])


AC_DEFUN([AC_EOSPAC5_FINALIZE], [dnl

   # set up the libraries
   if test -n "${vendor_eospac}"; then

       # set up library paths
       if test -z "${EOSPAC5_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_eospac, -l${with_eospac})
       elif test -n "${EOSPAC5_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_eospac, -L${EOSPAC5_LIB} -l${with_eospac})
       fi

       # add EOSPAC5 directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${EOSPAC5_LIB}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_LAPACK_SETUP
dnl
dnl LAPACK SETUP (on by default)
dnl LAPACK is a required vendor
dnl
dnl NOTE: this also sets up the BLAS
dnl
dnl note that we add system specific libraries to this list in
dnl ac_platforms.m4
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_LAPACK_SETUP], [dnl

   dnl define --with-lapack
   AC_ARG_WITH(lapack,
      [  --with-lapack=[vendor,atlas]
                          determine LAPACK implementation (vendor default)])

   dnl define --with-lapack-lib
   AC_WITH_DIR(lapack-lib, LAPACK_LIB, \${LAPACK_LIB_DIR}, 
	       [tell where LAPACK libs are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_lapack=$1

   # lapack is set to vendor by default
   if test "${with_lapack:=vendor}" = yes ; then
       with_lapack='vendor'
   fi

   # define the atlas libraries (these are system independent)
   if test "${with_lapack}" = atlas; then
       lapack_libs='-llapack -lf77blas -lcblas -latlas'
   fi
])


AC_DEFUN([AC_LAPACK_FINALIZE], [dnl

   # set up lapack libraries
   if test -n "${vendor_lapack}"; then

       # set libraries
       if test -z "${LAPACK_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_lapack, ${lapack_libs})
       elif test -n "${LAPACK_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_lapack, -L${LAPACK_LIB} ${lapack_libs})
       fi

       # add LAPACK directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${LAPACK_LIB}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_GRACE_SETUP
dnl
dnl GRACE SETUP (on by default)
dnl GRACE is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_GRACE_SETUP], [dnl

   dnl define --with-grace
   AC_ARG_WITH(grace,
      [  --with-grace=[lib]      determine the grace lib (grace_np is the default])
 
   dnl define --with-grace-inc
   AC_WITH_DIR(grace-inc, GRACE_INC, \${GRACE_INC_DIR},
	       [tell where GRACE includes are])

   dnl define --with-grace-lib
   AC_WITH_DIR(grace-lib, GRACE_LIB, \${GRACE_LIB_DIR},
	       [tell where GRACE libraries are])

   # set default value of grace includes and libs
   if test "${with_grace:=grace_np}" = yes ; then
       with_grace='grace_np'
   fi

   # define GRACE header file
   GRACE_H="<${with_grace}.h>"
   AC_DEFINE_UNQUOTED(GRACE_H, ${GRACE_H})dnl

   # determine if this package is needed for testing or for the 
   # package
   vendor_grace=$1

])


AC_DEFUN([AC_GRACE_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_grace}" ; then

       # include path
       if test -n "${GRACE_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${GRACE_INC}"
       fi

       # library path
       if test -n "${GRACE_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_grace, -L${GRACE_LIB} -l${with_grace})
       elif test -z "${GRACE_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_grace, -l${with_grace})
       fi

       # add GRACE directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${GRACE_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${GRACE_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_SPICA_SETUP
dnl
dnl SPICA LIBRARY SETUP (on by default -lSpicaCSG)
dnl SPICA is an optional vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_SPICA_SETUP], [dnl

   dnl define --with-spica
   AC_ARG_WITH(spica,
      [  --with-spica[=yes]                 spica is on by default])
	
   dnl define --with-spica-inc and --with-spica-lib
   AC_WITH_DIR(spica-inc, SPICA_INC, \${SPICA_INC_DIR},
	       [tell where SPICA includes are])
   AC_WITH_DIR(spica-lib, SPICA_LIB, \${SPICA_LIB_DIR},
	       [tell where SPICA libraries are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_spica=$1

   # define variable if spica is on
   if test "${with_spica:=yes}" != no; then
       AC_DEFINE([USE_SPICA])
   fi
])


AC_DEFUN([AC_SPICA_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_spica}"; then

       # include path
       if test -n "${SPICA_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${SPICA_INC}"
       fi
   
       # libraries
       if test -n "${SPICA_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_spica, -L${SPICA_LIB} -lSpicaCSG)
       elif test -z "${SPICA_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_spica, -lSpicaCSG)
       fi

       # add spica directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${SPICA_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${SPICA_INC}"

   fi
])

dnl-------------------------------------------------------------------------dnl
dnl AC_XERCES_SETUP
dnl
dnl XERCES LIBRARY SETUP
dnl xerces is a required vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_XERCES_SETUP], [dnl

   dnl define --with-xerces
   AC_ARG_WITH(xerces,
      [  --with-xerces[=lib]      determine the XERCES xml lib (xerces-c is default)])
	
   dnl define --with-xerces-inc and --with-xerces-lib
   AC_WITH_DIR(xerces-inc, XERCES_INC, \${XERCES_INC_DIR},
	       [tell where XERCES includes are])
   AC_WITH_DIR(xerces-lib, XERCES_LIB, \${XERCES_LIB_DIR},
	       [tell where XERCES libraries are])

   # determine if this package is needed for testing or for the 
   # package
   vendor_xerces=$1

   # default (xerces is set to xerces-c by default)
   if test "${with_xerces:=xerces-c}" = yes ; then
       with_xerces='xerces-c'
   fi
])


AC_DEFUN([AC_XERCES_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_xerces}"; then

       # include path
       if test -n "${XERCES_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${XERCES_INC}"
       fi
   
       # libraries
       if test -n "${XERCES_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_xerces, -L${XERCES_LIB} -l${with_xerces})
       elif test -z "${XERCES_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_xerces, -l${with_xerces})
       fi

       # add xerces directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${XERCES_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${XERCES_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_HDF5_SETUP
dnl
dnl HDF5 SETUP (on by default; 'mpi' if mpi in use, else 'serial')
dnl HDF5 is an optional vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_HDF5_SETUP], [dnl

   dnl define --with-hdf5
   AC_ARG_WITH(hdf5,
      [  --with-hdf5=[serial,mpi]      determine hdf5 implementation (default:  'mpi' if mpi in use, else 'serial')])
 
   dnl define --with-hdf5-inc
   AC_WITH_DIR(hdf5-inc, HDF5_INC, \${HDF5_INC_DIR},
	       [tell where HDF5 includes are])

   dnl define --with-hdf5-lib
   AC_WITH_DIR(hdf5-lib, HDF5_LIB, \${HDF5_LIB_DIR},
	       [tell where HDF5 libraries are])

   # default (mpi if mpi is in use, else serial)
   if test "${with_hdf5:=no}" = yes ; then
       if test "${with_mpi}" != no ; then
	   with_hdf5='mpi'
       else
	   with_hdf5='serial'
       fi
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_hdf5=$1

   # define variable if hdf5 is on
   if test "${with_hdf5:=yes}" != no; then
       AC_DEFINE([USE_HDF5])
   fi

])


AC_DEFUN([AC_HDF5_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_hdf5}" ; then

       # include path
       if test -n "${HDF5_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${HDF5_INC}"
       fi

       # library path
       if test -n "${HDF5_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_hdf5, -L${HDF5_LIB} -lhdf5 -lz)
       elif test -z "${HDF5_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_hdf5, -lhdf5 -lz)
       fi

       # add HDF5 directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${HDF5_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${HDF5_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_UDM_SETUP
dnl
dnl UDM SETUP (on by default; 'mpi' if mpi in use, else 'serial')
dnl UDM is an optional vendor
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_UDM_SETUP], [dnl

   dnl define --with-udm
   AC_ARG_WITH(udm,
      [  --with-udm=[serial,mpi]      determine udm implementation (default:  'mpi' if mpi in use, else 'serial')])
 
   dnl define --with-udm-inc
   AC_WITH_DIR(udm-inc, UDM_INC, \${UDM_INC_DIR},
	       [tell where UDM includes are])

   dnl define --with-udm-lib
   AC_WITH_DIR(udm-lib, UDM_LIB, \${UDM_LIB_DIR},
	       [tell where UDM libraries are])

   # default (mpi if mpi is in use, else serial)
   if test "${with_udm:=no}" = yes ; then
       if test "${with_mpi}" != no ; then
	   with_udm='mpi'
       else
	   with_udm='serial'
       fi
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_udm=$1

   # define variable if udm is on
   if test "${with_udm:=no}" != no; then
       AC_DEFINE([USE_UDM])
   fi

])


AC_DEFUN([AC_UDM_FINALIZE], [dnl

   # set up the libraries and include path
   if test -n "${vendor_udm}" ; then

       # include path
       if test -n "${UDM_INC}"; then
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${UDM_INC}"
           # set extra #define if using udm in parallel
           if test "${with_udm}" = mpi ; then
               AC_DEFINE(UDM_HAVE_PARALLEL)
           fi
       fi

       # library path
       if test -n "${UDM_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_udm, -L${UDM_LIB} -ludm)
       elif test -z "${UDM_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_udm, -ludm)
       fi

       # add UDM directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${UDM_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${UDM_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_DLOPEN_SETUP
dnl
dnl This is an optional vendor.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DLOPEN_SETUP], [dnl

   dnl define --enable-dlopen
   AC_ARG_ENABLE(dlopen,
      [  --enable-dlopen          Enable dlopen (default: on if --enable-shared, off otherwise)])

   # determine if this package is needed for testing or for the
   # package.
   vendor_dlopen=$1 

   # set default value for enable_dlopen, which is the value of enable_shared.
   if test "${enable_shared}" = yes ; then
       if test "${enable_dlopen:=yes}" != no ; then 
	   enable_dlopen=yes
       fi
   else
       if test "${enable_dlopen:=no}" != no ; then 
	   enable_dlopen=yes
       fi
   fi

   # turn off dlopen if not using shared libraries.
   if test "${enable_shared}" != yes ; then
       if test "${enable_dlopen}" = yes ; then
	   AC_MSG_WARN("Must specify --enable-shared when using --enable-dlopen.")
           AC_MSG_WARN("   dlopen disabled.")
       fi
       enable_dlopen=no
   fi

   if test "${enable_dlopen}" = yes ; then
       AC_DEFINE(USE_DLOPEN)
   fi
]) 


AC_DEFUN([AC_DLOPEN_FINALIZE], [dnl
   # Libraries are platform-specific; done in ac_platforms.
])

dnl-------------------------------------------------------------------------dnl
dnl AC_VENDOR_FINALIZE
dnl
dnl Run at the end of the environment setup to add defines required by
dnl the vendors.  We do this to allow platform specific mods to the 
dnl vendor defines BEFORE they are added to CCPFLAGS, etc. 
dnl
dnl This macro needs to be updated when new vendors are added.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_VENDOR_FINALIZE], [dnl

   # call finalize functions for each vendor, the order is important
   # each vendor setup is appended to the previous; thus, the calling
   # level goes from high to low
   AC_TRILINOS_FINALIZE
   AC_GSL_FINALIZE

   AC_AZTEC_FINALIZE
   AC_PCG_FINALIZE
   AC_HYPRE_FINALIZE
   AC_SCALAPACK_FINALIZE
   AC_BLACS_FINALIZE
   AC_LAPACK_FINALIZE
   AC_EOSPAC5_FINALIZE
   AC_GANDOLF_FINALIZE
   AC_SPRNG_FINALIZE
   AC_GRACE_FINALIZE
   AC_METIS_FINALIZE
   AC_SPICA_FINALIZE
   AC_XERCES_FINALIZE

   AC_UDM_FINALIZE
   AC_HDF5_FINALIZE

   AC_MPI_FINALIZE
   AC_DLOPEN_FINALIZE

   # print out vendor include paths
   AC_MSG_CHECKING("vendor include paths")
   if test -n "${VENDOR_INC_DIRS}"; then
       AC_MSG_RESULT("${VENDOR_INC_DIRS}")
   else
       AC_MSG_RESULT("no vendor include dirs defined")
   fi

   # print out vendor lib paths
   AC_MSG_CHECKING("vendor lib paths")
   if test -n "${VENDOR_LIB_DIRS}"; then
       AC_MSG_RESULT("${VENDOR_LIB_DIRS}")
   else
       AC_MSG_RESULT("no vendor lib dirs defined")
   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_ALL_VENDORS_SETUP
dnl
dnl DRACO INCLUSIVE VENDOR MACRO
dnl-------------------------------------------------------------------------dnl
dnl allows one to include all vendor macros by calling this macro.
dnl designed for draco/configure.in and draco/src/configure.in

AC_DEFUN(AC_ALL_VENDORS_SETUP, [dnl

   dnl include all macros for easy use in top-level configure.in's
   AC_MPI_SETUP(pkg)
   AC_SPRNG_SETUP(pkg)
   AC_PCG_SETUP(pkg)
   AC_AZTEC_SETUP(pkg)
   AC_GSL_SETUP(pkg)
   AC_TRILINOS_SETUP(pkg)
   AC_METIS_SETUP(pkg)
   AC_LAPACK_SETUP(pkg)
   AC_GANDOLF_SETUP(pkg)
   AC_EOSPAC5_SETUP(pkg)
   AC_GRACE_SETUP(pkg)
   AC_SPICA_SETUP(pkg)
   AC_XERCES_SETUP(pkg)
   AC_HDF5_SETUP(pkg)
   AC_UDM_SETUP(pkg)
   AC_DLOPEN_SETUP(pkg)
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_vendors.m4
dnl-------------------------------------------------------------------------dnl


dnl-------------------------------------------------------------------------dnl
dnl File  : draco/config ac_dracoenv.m4
dnl Author: Thomas M. Evans
dnl Date  : 1999/02/04 01:56:21
dnl
dnl Defines the Draco build system environment.  This is the main
dnl configure function.
dnl
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_STLPORT_ENV
dnl
dnl Used by AC_DRACO_ENV, this macro checks the configure line for the
dnl presence of "--with-stlport".  If this option is found, the build
dnl system's environment is modified so that all the all C++ compiles
dnl use the STL libraries included with STLPort instead of the
dnl compiler's native STL defintions.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DBS_STLPORT_ENV, [dnl

   AC_MSG_CHECKING("for stlport")
   if test "${with_stlport:=no}" != no; then
      if ! test -d "${with_stlport}/include/stlport"; then
         AC_MSG_ERROR("Invalid directory ${with_stlport}/include/stlport")
      fi
      CPPFLAGS="-I${with_stlport}/include/stlport ${CPPFLAGS}"
      CXXFLAGS="-I${with_stlport}/include/stlport ${CXXFLAGS}"
      AC_MSG_RESULT("-I${with_stlport}/include added to CXXFLAGS.")
      case $with_cxx in
      sgi)
         stlport_libname='mipspro'
         stlport_xlinker=' '
         ;;
      gcc) dnl for everything else use gcc
         stlport_libname='gcc'
         stlport_xlinker="-Xlinker"
         ;;
      *) 
         AC_MSG_ERROR("stlport not available with this compiler.")
         ;;
      esac
      AC_MSG_CHECKING("for debug stlport mode")
      if test "${enable_debug:-yes}" = yes; then
         if ! test -r "${with_stlport}/lib/libstlport_${stlport_libname}_stldebug.a"; then
            AC_MSG_ERROR("Invalid library ${with_stlport}/lib/libstlport_${stlport_libname}_stldebug.a")
         fi
         LIBS="-L${with_stlport}/lib -lstlport_${stlport_libname}_stldebug ${LIBS}"
         CXXFLAGS="${CXXFLAGS} -D_STLP_DEBUG"
         AC_MSG_RESULT("yes")
      else
         if ! test -r "${with_stlport}/lib/libstlport_${stlport_libname}.a"; then
            AC_MSG_ERROR("Invalid library ${with_stlport}/lib/libstlport_${stlport_libname}.a")
         fi
         LIBS="-L${with_stlport}/lib -lstlport_${stlport_libname} ${LIBS}"
         AC_MSG_RESULT("no")
      fi

      AC_MSG_CHECKING("for RPATH mods for stlport")
      RPATH="${stlport_xlinker} -rpath ${with_stlport}/lib ${RPATH}"
      AC_MSG_RESULT("Added ${stlport_xlinker} -rpath ${with_stlport}/lib to RPATH")

   elif test "${with_stlport}" = yes; then
      AC_MSG_ERROR("Must define path to stlport when using --with-stlport=[dir]")
   else
      AC_MSG_RESULT("none")
   fi

   dnl end of AC_DBS_STLPORT_ENV
])

dnl-------------------------------------------------------------------------dnl
dnl AC_DRACO_ENV
dnl
dnl Assembles the Draco build system compile-time environment.  
dnl It processes all of the options given to configure.  It does
dnl NOT do any compile or link testing.  That functionality is
dnl defined in ac_dracotests.m4.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_ENV, [dnl

   dnl
   dnl CONFIGURE ARGUMENTS
   dnl

   # Retrieve the configure command line for possible use in 
   # regression test output.

   configure_command="$[]0 $[]*"

   dnl
   dnl ADD DRACO CONFIGURE ARGUMENTS
   dnl

   AC_DRACO_ARGS

   dnl
   dnl first find the host
   dnl
   
   AC_REQUIRE([AC_CANONICAL_HOST])

   dnl
   dnl INSTALL
   dnl

   # we use the install script provided with autoconf on all machines
   INSTALL='${config_dir}/install-sh -c'
   INSTALL_DATA='${INSTALL} -m 644'

   dnl
   dnl C4 OPERATIONS
   dnl

   # do the correct #defines
   if test "$with_c4" = scalar ; then
       AC_DEFINE(C4_SCALAR)
   elif test "$with_c4" = mpi ; then
       AC_DEFINE(C4_MPI)
   fi

   dnl
   dnl DBC SETUP
   dnl

   # set the DBC level
   if test "${with_dbc:=default}" != default ; then
       AC_DEFINE_UNQUOTED(DBC, $with_dbc)
   fi

   dnl
   dnl LIBRARIES
   dnl
   
   # set the libsuffix variable
   if test "${enable_shared:=no}" = yes ; then
       libsuffix='.so'
   else
       libsuffix='.a'
   fi

   dnl      
   dnl POSIX SOURCE
   dnl

   dnl system dependent posix defines are performed in the
   dnl SYSTEM-SPECIFIC SETUP section below

   dnl
   dnl TOOL CHECKS 
   dnl

   # the tool checks are called in the top-level configure, so in 
   # each subsequent configure these should just grab cached values
   AC_DRACO_CHECK_TOOLS dnl

   dnl
   dnl COMPILER SETUPS
   dnl

   # the default compiler is C++; we do not turn on F90 unless
   # AC_WITH_F90 is called in configure.in (which sets with_cxx='no')
   if test "${with_cxx}" = no ; then

       # if with_f90 defined test with_f90 for compiler, and call setup
       # if with_f90 set to yes or not set 
       # attempt to guess compiler based on target
       AC_F90_ENV dnl

   else
   
       # set up the C++ compilers; if with_cxx is undefined, an
       # appropriate default for the machine will be choosen
       AC_CPP_ENV dnl

   fi

   # STL port checks and setup

   AC_DBS_STLPORT_ENV

   dnl add any additional flags

   # add user defined cppflags
   if test "${with_cppflags:=no}" != no ; then
       CPPFLAGS="${with_cppflags} ${CPPFLAGS}"
   fi

   # add user defined cxxflags
   if test "${with_cxxflags:=no}" != no ; then
       CXXFLAGS="${with_cxxflags} ${CXXFLAGS}"
   fi

   # add user defined cflags
   if test "${with_cflags:=no}" != no ; then
       CFLAGS="${with_cflags} ${CFLAGS}"
   fi

   # add user defined f90flags
   if test "${with_f90flags:=no}" != no ; then
       F90FLAGS="${with_f90flags} ${F90FLAGS}"
   fi

   # add user defined ARFLAGS
   if test "${with_arflags:=no}" != no ; then
       ARFLAGS="${with_arflags} ${ARFLAGS}"
   fi

   # add user defined LDFLAGS
   if test "${with_ldflags:=no}" != no ; then
       LDFLAGS="${with_ldflags} ${LDFLAGS}"
   fi

   # check user added libs (using --with-libs); these are appended to
   # LIBS after the machine-specific setup
   if test "${with_libs}" = yes ; then
       AC_MSG_ERROR("Must define libs when using --with-libs")
   fi

   dnl throw message errors for poorly defined flags
   
   if test "${with_cxxflags}" = yes || test "${with_cflags}" = yes ||\
      test "${with_f90flags}" = yes || test "${with_arflags}" = yes \
      || test "${with_ldflags}" = yes \
      || test "${with_cppflags}" = yes ; then
       AC_MSG_ERROR("Poor definition of user defined flags!")
   fi
   
   dnl check for ranlib
   AC_PROG_RANLIB

   dnl
   dnl SYSTEM-SPECIFIC SETUP
   dnl

   # this function macro sets up all of the platform specific 
   # environment parameters (except compilers)
   AC_DBS_PLATFORM_ENVIRONMENT dnl

   # add user-defined libraries
   LIBS="${LIBS} ${with_libs} -lm"

   dnl
   dnl DRACO TEST SYSTEM
   dnl

   # determine whether this is a scalar or parallel test suite,
   # the tests can be inherently scalar or they can be the result
   # of a parallel build

   test_scalar='no'

   # If we ran AC_RUNTESTS with "serial" then mark it so here.
   for np in $test_nprocs; do
       if test $np = serial || test $np = scalar ; then
          test_scalar="scalar"
       fi
   done

   # if this is a parallel build, mark the tests scalar
   if test "${with_c4}" = scalar ; then
       test_scalar="scalar"
   fi

   # define the TESTFLAGS, for parallel runs the processor will be
   # added later in the Makefile

   if test "${test_scalar}" = scalar ; then
       test_flags="--${test_exe:=binary}"
   elif test "${with_c4}" = mpi ; then
       test_flags="--${test_exe:=binary} --mpi"
   fi

   ## define the test_output_files for cleaning
   for file in $test_alltarget; do
       if test "${test_scalar}" = scalar ; then
	   test_output_files="${test_output_files} ${file}-scalar.log"
       else
	   for np in $test_nprocs; do
	       test_output_files="${test_output_files} ${file}-${np}.log"
	   done
       fi
   done

   # Define the package-level source directory (e.g. draco)
   AC_FIND_TOP_SRC($srcdir, package_top_srcdir)

   dnl
   dnl ENVIRONMENT SUBSTITUTIONS
   dnl

   AC_DBS_VAR_SUBSTITUTIONS

   dnl end of AC_DRACO_ENV
])


dnl-------------------------------------------------------------------------dnl
dnl end of ac_dracoenv.m4
dnl-------------------------------------------------------------------------dnl


dnl-------------------------------------------------------------------------dnl
dnl ac_dracoarg.m4
dnl
dnl Declarations of Draco configure options (with some default
dnl settings). 
dnl
dnl Thomas M. Evans
dnl 1999/02/04 01:56:20
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl AC_DRACO_ARGS
dnl
dnl Declaration of Draco non-vendor configure options. This macro can 
dnl be called to fill out configure help screens
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_ARGS, [dnl

   dnl
   dnl Library prefix
   dnl
     
   AC_ARG_WITH(lib-prefix,
      [  --with-lib-prefix[=library prefix]
                          give prefix to libraries (default rtt_)])

   # default for lib_prefix is rtt_
   LIB_PREFIX="${with_lib_prefix:=rtt_}"
   if test "${LIB_PREFIX}" = no ; then
       LIB_PREFIX=''
   fi

   dnl
   dnl c4 toggle (scalar by default)
   dnl

   dnl define --with-c4
   AC_ARG_WITH(c4, 
      [  --with-c4[=scalar,mpi,shmem]   
		          turn on c4 (default scalar) ])

   # give with-c4 implied argument
   if test "${with_c4:=scalar}" = yes ; then
       with_c4='scalar'
   fi

   dnl
   dnl DBC toggle
   dnl

   dnl defines --with-dbc
   AC_ARG_WITH(dbc,
      [  --with-dbc[=level]      set Design-by-Contract])
	
   if test "${with_dbc}" = yes ; then
       with_dbc='7'
   elif test "${with_dbc}" = no ; then
       with_dbc='0'
   fi
	
   dnl
   dnl SHARED versus ARCHIVE libraries
   dnl

   dnl defines --enable-shared
   AC_ARG_ENABLE(shared,
      [  --enable-shared         turn on shared libraries (.a default)])

   dnl
   dnl CHOOSE A C++ COMPILER
   dnl

   dnl defines --with-cxx
   AC_ARG_WITH(cxx,
      [  --with-cxx[=gcc,icpc,sgi,kcc,compaq,guide]                                    
                          choose a c++ compiler (defaults are machine dependent)])

   dnl the default is gcc
   if test "${with_cxx}" = yes ; then
       with_cxx='gcc'
   fi

   dnl
   dnl STATIC VERSUS DYNAMIC LINKING
   dnl

   dnl defines --enable-static-ld
   AC_ARG_ENABLE(static-ld,
      [  --enable-static-ld      use (.a) libraries if possible])

   dnl
   dnl ANSI STRICT COMPLIANCE
   dnl

   dnl defines --enable-strict-ansi
   AC_ARG_ENABLE(strict-ansi,
      [  --disable-strict-ansi   turn off strict ansi compliance])

   dnl
   dnl ONE_PER INSTANTIATION FLAG
   dnl

   dnl defines --enable-one-per
   AC_ARG_ENABLE(one-per,
      [  --disable-one-per       turn off --one_per flag])

   dnl
   dnl COMPILER OPTIMZATION LEVEL
   dnl

   dnl defines --with-opt
   AC_ARG_WITH(opt,
      [  --with-opt[=0,1,2,3]    set optimization level (0 by default)])

   if test "${with_opt}" = yes ; then
       with_opt='0'
   fi

   dnl defines --enable-debug
   AC_ARG_ENABLE(debug,
      [  --enable-debug          turn on debug (-g) option])

   dnl
   dnl POSIX SOURCE
   dnl

   dnl defines --with-posix
   AC_ARG_WITH(posix,
      [  --with-posix[=num]      give posix source (system-dependent defaults)])

   dnl
   dnl ADD TO CPPFLAGS
   dnl
   
   dnl defines --with-cppflags
   AC_ARG_WITH(cppflags,
      [  --with-cppflags[=flags] add flags to \$CPPFLAGS])

   dnl
   dnl ADD TO CXXFLAGS
   dnl
   
   dnl defines --with-cxxflags
   AC_ARG_WITH(cxxflags,
      [  --with-cxxflags[=flags] add flags to \$CXXFLAGS])

   dnl
   dnl ADD TO CFLAGS
   dnl
   
   dnl defines --with-cflags
   AC_ARG_WITH(cflags,
      [  --with-cflags[=flags]   add flags to \$CFLAGS])

   dnl
   dnl ADD TO F90FLAGS
   dnl
   
   dnl defines --with-f90flags
   AC_ARG_WITH(f90flags,
      [  --with-f90flags[=flags] add flags to \$F90FLAGS])

   dnl
   dnl ADD TO ARFLAGS
   dnl
   
   dnl defines --with-arflags
   AC_ARG_WITH(arflags,
      [  --with-arflags[=flags]  add flags to \$ARFLAGS])

   dnl
   dnl ADD TO LDFLAGS
   dnl
   
   dnl defines --with-ldflags
   AC_ARG_WITH(ldflags,
      [  --with-ldflags[=flags]  add flags to \$LDFLAGS])

   dnl 
   dnl ADD TO LIBRARIES
   dnl

   dnl defines --with-libs
   AC_ARG_WITH(libs,
      [  --with-libs=[libs]      add libs to \$LIBS])

   dnl
   dnl CHOSE BIT COMPILATION ON SGI'S
   dnl

   dnl defines --enable-32-bit
   AC_ARG_ENABLE(32-bit,
      [  --enable-32-bit         do 32-bit compilation (compiler dependent)])

   dnl defines --enable-64-bit
   AC_ARG_ENABLE(64-bit,
      [  --enable-64-bit         do 64-bit compilation (compiler dependent)])

   dnl
   dnl CHOSE MIPS INSTRUCTION SET ON SGI'S
   dnl

   dnl defines --with-mips
   AC_ARG_WITH(mips,
      [  --with-mips[=1,2,3,4]   set mips, mips4 by default (SGI ONLY)])

   if test "${with_mips}" = yes ; then
       with_mips='4'
   fi

   dnl 
   dnl STLport
   dnl

   dnl specify location of stlport installation.
   AC_ARG_WITH(stlport,
      [  --with-stlport        replace default STL with stlPort (off by default)])

   dnl Doxygen options

   AC_ARG_ENABLE(latex-doc,
      [  --enable-latex-doc      build latex docs with doxygen (off by default)],
      [AC_SUBST(latex_yes_no,'YES')],
      [AC_SUBST(latex_yes_no,'NO')])

   AC_ARG_WITH(doc-output,
      [  --with-doc-output=path  build documentation in path (prefix/documentation by default)],
      [AC_SUBST(doxygen_output_top,${with_doc_output})],
      [doxygen_output_top='DEFAULT'])

   dnl end of AC_DRACO_ARGS
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_dracoarg.m4
dnl-------------------------------------------------------------------------dnl


dnl ========================================================================
dnl 
dnl 	Author:	Mark G. Gray
dnl 		Los Alamos National Laboratory
dnl 	Date:	Wed Apr 19 16:39:19 MDT 2000
dnl 
dnl 	Copyright (c) 2000 U. S. Department of Energy. All rights reserved.
dnl 
dnl ========================================================================

dnl NAME

dnl	AC_WITH_F90, AC_F90_ENV

dnl SYNOPSIS/USAGE

dnl     AC_WITH_F90
dnl     AC_F90_ENV

dnl DESCRIPTION

dnl     AC_WITH_F90 sets the variable with_f90 to yes if it is not already 
dnl     set.

dnl     AC_F90_ENV set environment variables F90, F90FLAGS, F90EXT, 
dnl     F90FREE, F90FIXED, and MODFLAG for the compiler requested by 
dnl     with_f90.  If no specific compiler is requested, guess a compiler 
dnl     based on the target
dnl
========================================================================

dnl ### Ensure with_f90 set
AC_DEFUN(AC_WITH_F90, [dnl
   : ${with_f90:=yes}
    
   dnl turn off C++ compiler
   with_cxx='no'

   dnl defines --with-f90
   AC_ARG_WITH(f90,
       [  --with-f90[=XL,Fujitsu,Lahey,Portland,WorkShop,Cray,MIPS,Compaq,HP,Intel,NAG,Absoft]
                          choose an F90 compiler])
])

dnl
dnl CHOOSE A F90 COMPILER
dnl

AC_DEFUN(AC_F90_ENV, [dnl
   AC_REQUIRE([AC_CANONICAL_HOST])

   case "${with_f90:=yes}" in
   XL)
       AC_COMPILER_XL_F90
   ;;
   Fujitsu)
       AC_COMPILER_FUJITSU_F90
   ;;
   Lahey)
       AC_COMPILER_LAHEY_F90
   ;;
   Portland)
       AC_COMPILER_PORTLAND_F90
   ;;
   WorkShop)
       AC_COMPILER_WORKSHOP_F90
   ;;
   Cray)
      AC_COMPILER_CRAY_F90
   ;;
   MIPS)
       AC_COMPILER_MIPS_F90
   ;;
   Compaq)
       AC_COMPILER_COMPAQ_F90
   ;;
   HP)
       AC_COMPILER_HP_F90
   ;;
   Intel)
       AC_COMPILER_INTEL_F90
   ;;
   NAG)
       AC_COMPILER_NAG_F90
   ;;
   Absoft)
       AC_COMPILER_ABSOFT_F90
   ;;
   yes)				# guess compiler from target platform
       case "${host}" in   
       rs6000-ibm-aix*)
           AC_COMPILER_XL_F90
       ;;
       powerpc-ibm-aix*)
           AC_COMPILER_XL_F90
       ;;
       sparc-sun-solaris2.*)
           AC_COMPILER_WORKSHOP_F90
       ;;
       i?86-pc-linux*)
           AC_COMPILER_LAHEY_F90
       ;;
       ymp-cray-unicos*)
          AC_COMPILER_CRAY_F90
       ;;
       mips-sgi-irix*)
          AC_COMPILER_MIPS_F90
       ;;
       i??86-pc-cygwin*)
          AC_COMPILER_COMPAQ_F90
       ;;
       alpha*)
          AC_COMPILER_COMPAQ_F90
       ;;
       *hp-hpux*)
          AC_COMPILER_HP_F90
       ;;
       *)
          AC_MSG_ERROR([Cannot guess F90 compiler, set --with-f90])
       ;;
       esac
   ;;
   no)
   ;;
   *)
       AC_MSG_ERROR([Unrecognized F90 compiler, use --help])
   ;;
   esac

   AC_SUBST(F90FREE)
   AC_SUBST(F90FIXED)
   AC_SUBST(F90FLAGS)
   AC_SUBST(MODFLAG)
])

dnl-------------------------------------------------------------------------dnl
dnl IBM XLF95 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_XL_F90, [dnl

   # Check for working XL F90 compiler

  if test "${with_upslib:=no}" != "no"
  then
     AC_CHECK_PROG(F90, mpxlf95, mpxlf95, none)
     if test "${F90}" != mpxlf95
     then
         AC_MSG_ERROR([not found])
     fi
  else
     AC_CHECK_PROG(F90, xlf95, xlf95, none)
     if test "${F90}" != xlf95
     then
         AC_MSG_ERROR([not found])
     fi
  fi
  
   # FREE, FIXED AND MODULE FLAGS

   F90FREE='-qfree=f90'
   F90FIXED='-qfixed'
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=

   # COMPILATION FLAGS

   if test "$F90FLAGS" = ""
   then
     # F90FLAGS="-qsuffix=f=f90 -qmaxmem=-1 -qextchk -qarch=pwr2 -bmaxstack:0x70000000 -bmaxdata:0x70000000 -qalias=noaryovrlp -qhalt=s ${F90FREE}"
       F90FLAGS="-qsuffix=f=f90 -qmaxmem=-1 -qextchk -qarch=auto -bmaxstack:0x70000000 -bmaxdata:0x70000000 -qalias=noaryovrlp -qnosave -qlanglvl=95pure -qzerosize ${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	   trapflags="-qinitauto=FF"
	   trapflags="${trapflags} -qflttrap=overflow:underflow:zerodivide:invalid:enable"
	   trapflags="${trapflags} -qsigtrap"
	   F90FLAGS="-g -d -C ${trapflags} -bloadmap:loadmap.dat ${F90FLAGS}"
       else
	 # F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
	   F90FLAGS="-O3 ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_XL_F90
])

dnl-------------------------------------------------------------------------dnl
dnl FUJITSU F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_FUJITSU_F90, [dnl

   # Check for working Fujitsu F90 compiler

   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90 && ${F90} -V 2>&1 | grep "Fujitsu"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG

   F90FREE='-Free'
   F90FIXED='-Fixed'
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-static-flib'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="-X9 -Am ${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g -Haesu ${F90FLAGS}"
       else
	    F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_FUJITSU_F90
])

dnl-------------------------------------------------------------------------dnl
dnl LAHEY F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_LAHEY_F90, [dnl

   AC_CHECK_PROG(F90, lf95, lf95, none)

   # F90FREE, F90FIXED AND MODFLAG

   F90FREE='--nfix'
   F90FIXED='--fix'
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-static-flib'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
     # F90FLAGS="--f95 ${F90FREE}"
       F90FLAGS="--staticlink --f95 --in --info --swm 2004,2006,2008,8202,8203,8204,8205,8206,8209,8220 ${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	  # F90FLAGS="-g --chk --trace ${F90FLAGS}"
	    F90FLAGS="-g --ap --chk --pca --private --trap --wo ${F90FLAGS}"
       else
	  # F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
	    F90FLAGS="-O --ntrace ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_LAHEY_F90
])

dnl-------------------------------------------------------------------------dnl
dnl PORTLAND F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_PORTLAND_F90, [dnl

   # Check for working Portland Group F90 compiler

   AC_CHECK_PROG(F90, pgf90, pgf90, none)
   if test "${F90}" = pgf90 && ${F90} --V 2>&1 | grep "Portland"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG

   F90FREE='-Mfreeform'
   F90FIXED='-Mnofreeform'
   MODFLAG='-module'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC=

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g -Mbounds -Mchkptr ${F90FLAGS}"
       else
	    F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_PORTLAND_F90
])

dnl-------------------------------------------------------------------------dnl
dnl COMPAQ F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_COMPAQ_F90, [dnl

   # Check for working compaq F90 compiler

   AC_CHECK_PROG(F90, f95, f95, none)
   if test "${F90}" = f95 && ${F90} -version 2>&1 | grep "Fortran"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG

   F90FREE=''
   F90FIXED=''
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-non_shared'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
     # F90FLAGS="${F90FREE} -assume byterecl"
       F90FLAGS="${F90FREE} -assume byterecl -automatic -std -warn argument_checking"

       if test "${enable_debug:=no}" = yes
       then
	  # F90FLAGS="-g ${F90FLAGS}"
	    F90FLAGS="-g -check bounds -fpe2 ${F90FLAGS}"
       else
	  # F90FLAGS="-O ${F90FLAGS}"
	    F90FLAGS="-O5 -arch host -assume noaccuracy_sensitive -math_library accurate -tune host ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_COMPAQ_F90
])

dnl-------------------------------------------------------------------------dnl
dnl SUN WORKSHOP F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_WORKSHOP_F90, [dnl

   # Check for working WorkShop F90 compiler

   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90 && ${F90} -V 2>&1 | grep "WorkShop"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # Set F90FREE, F90FIXED, and MODFLAG

   F90FREE='-free'
   F90FIXED='-fixed'
   MODFLAG='-M'

   # Set LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-Bstatic'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g"
       else
	    F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_WORKSHOP_F90
])

dnl-------------------------------------------------------------------------dnl
dnl CRAY_F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_CRAY_F90, [dnl

   # Check for working Cray F90 compiler

   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # FREE, FIXED AND MODULE FLAGS

   F90FREE='-f free'
   F90FIXED='-f fixed'
   MODFLAG='-p'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE}"

       if test "${enable_debug:=no}" = yes
       then
	   F90FLAGS="-g ${F90FLAGS}"
       else
	   F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_CRAY_F90
])

dnl-------------------------------------------------------------------------dnl
dnl IRIX MIPS F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_MIPS_F90, [dnl

   # Look for working MIPS compiler

   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90 && ${F90} -version 2>&1 | grep "MIPS"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # Set F90FREE, F90FIXED, and MODFLAG

   F90FREE='-freeform'
   F90FIXED='-col72'
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)

   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
	#F90FLAGS="${F90FREE} -OPT:Olimit=0"
	F90FLAGS="${F90FREE} -mips4 -r10000 -DEBUG:fullwarn=ON:woff=878,938,1193,1438"

	if test "${enable_debug:=no}" = yes
	then
	  # F90FLAGS="-g ${F90FLAGS}"
	    F90FLAGS="-g -check_bounds -DEBUG:trap_uninitialized=ON ${F90FLAGS}"
	else
	  # F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
	    F90FLAGS="-O3 -OPT:IEEE_arithmetic=2:roundoff=2 ${F90FLAGS}"
	fi
   fi

   dnl end of AC_COMPILER_MIPS_F90
])

dnl-------------------------------------------------------------------------dnl
dnl HP F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_HP_F90, [dnl

   # CHECK FOR WORKING HP F90 COMPILER
   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90 && ${F90} +version 2>&1 | grep "HP"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG
   F90FREE='+source=free'
   F90FIXED='+source=fixed'
   MODFLAG='-I'

   # LINKER AND LIBRARY (AR)
   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='+noshared'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE} +U77"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g -C ${F90FLAGS}"
       else
	    F90FLAGS="-O${with_opt:=} ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_HP_F90
])

dnl-------------------------------------------------------------------------dnl
dnl INTEL F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_INTEL_F90, [dnl

   # CHECK FOR WORKING INTEL F90 COMPILER
   AC_CHECK_PROG(F90, ifc, ifc, none)
   if test "${F90}" = ifc && ${F90} -V 2>&1 | grep "Intel"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG
   F90FREE='-FR'
   F90FIXED='-FI'
   MODFLAG='-I '
   MODSUFFIX='mod'

   # LINKER AND LIBRARY (AR)
   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-static'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE} -e95"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g -C -implicitnone ${F90FLAGS}"
       else
	    F90FLAGS="-O3 -fno-alias -tpp7 -ipo -pad -align ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_INTEL_F90
])

dnl-------------------------------------------------------------------------dnl
dnl NAG F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_NAG_F90, [dnl

   # CHECK FOR WORKING NAG F90 COMPILER
   AC_CHECK_PROG(F90, f95, f95, none)
   if test "${F90}" = f95 && ${F90} -V 2>&1 | grep "NAGWare"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
   # F90FREE, F90FIXED AND MODFLAG
   F90FREE='-free'
   F90FIXED='-fixed'
   MODFLAG='-I '

   # LINKER AND LIBRARY (AR)
   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC='-unsharedf95'

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="${F90FREE} -colour -info -target=native"

       if test "${enable_debug:=no}" = yes
       then
          # only use first line if memory error is suspected, too much output
          #   otherwise
	  # F90FLAGS="-g -C -mtrace=size -nan -u ${F90FLAGS}"
	    F90FLAGS="-g -C -nan -u ${F90FLAGS}"
       else
	    F90FLAGS="-O4 ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_NAG_F90
])

dnl-------------------------------------------------------------------------dnl
dnl ABSOFT F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_COMPILER_ABSOFT_F90, [dnl

   # CHECK FOR WORKING ABSOFT F90 COMPILER
   AC_CHECK_PROG(F90, f95, f95, none)
  
   # F90FREE, F90FIXED AND MODFLAG
   F90FREE=''
   F90FIXED=''
   MODFLAG='-p '

   # LINKER AND LIBRARY (AR)
   LD='${F90}'
   AR='ar'
   ARFLAGS=
   ARLIBS=
   F90STATIC=''

   # SET COMPILATION FLAGS IF NOT SET IN ENVIRONMENT
   if test "$F90FLAGS" = ""
   then
       F90FLAGS="-cpu:host -en"

       if test "${enable_debug:=no}" = yes
       then
	    F90FLAGS="-g -et -m0 -M399,1193,878 -Rb -Rc -Rs -Rp -trap=ALL ${F90FLAGS}"
       else
	    F90FLAGS="-O3 ${F90FLAGS}"
       fi
   fi

   dnl end of AC_COMPILER_ABSOFT_F90
])

dnl ========================================================================

dnl-------------------------------------------------------------------------dnl
dnl ac_compiler.m4
dnl
dnl Sets up all of the C++ compilers.
dnl
dnl Thomas M. Evans
dnl 1999/03/05 18:16:55
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl C++ COMPILER SETUP FUNCTION-->this is called within AC_DRACO_ENV;
dnl default is to use the C++ Compiler.  To change defaults,
dnl AC_WITH_F90 should be called in configure.in (before
dnl AC_DRACO_ENV)
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_CPP_ENV, [dnl

   # make sure that the host is defined
   AC_REQUIRE([AC_CANONICAL_HOST])

   dnl set up a default compiler
   case $host in

   # IRIX -> CC
   mips-sgi-irix6.*)   
       if test -z "${with_cxx}" ; then
	   with_cxx='sgi'
       fi
   ;;

   # COMPAQ -> CXX
   alpha*-dec-osf*)
       if test -z "${with_cxx}" ; then
	   with_cxx='compaq'
       fi
   ;;

   # IBM ASCI WHITE -> newxlC (use --with-cxx=ibm for regular SP2)
   *ibm-aix*)
       if test -z "${with_cxx}" ; then
	   with_cxx='asciwhite'
       fi
   ;;

   # EVERYTHING ELSE -> gcc
   *)
       if test -z "${with_cxx}" ; then
	   with_cxx='gcc'
       fi
   ;;
   esac

   dnl determine which compiler we are using

   # do tests of --with-cxx, see if the compiler exists and then call
   # the proper setup function
   
   if test "${with_cxx}" = sgi ; then
       AC_CHECK_PROG(CXX, CC, CC)
       AC_CHECK_PROG(CC, cc, cc)  

       if test "${CXX}" = CC && test "${CC}" = cc ; then
	   AC_DRACO_SGI_CC
       else 
	   AC_MSG_ERROR("Did not find SGI CC compiler!")
       fi

   elif test "${with_cxx}" = gcc ; then 
       AC_CHECK_PROG(CXX, g++, g++)
       AC_CHECK_PROG(CC, gcc, gcc)

       if test "${CXX}" = g++ && test "${CC}" = gcc ; then
	   AC_DRACO_GNU_GCC
       else
	   AC_MSG_ERROR("Did not find gnu c++ compiler!")
       fi

   elif test "${with_cxx}" = compaq ; then
       AC_CHECK_PROG(CXX, cxx, cxx)
       AC_CHECK_PROG(CC, cc, cc)

       if test "${CXX}" = cxx && test "${CC}" = cc ; then
	   AC_DRACO_COMPAQ_CXX
       else
	   AC_MSG_ERROR("Did not find Compaq cxx compiler!")
       fi

   elif test "${with_cxx}" = icpc ; then 
       AC_CHECK_PROG(CXX, icpc, icpc)

       if test "${CXX}" = icpc ; then
	   CC='icpc'
	   AC_DRACO_INTEL_ICPC
       else
	   AC_MSG_ERROR("Did not find Intel icpc compiler!")
       fi

   elif test "${with_cxx}" = pgi ; then
       # only allow PGI on LINUX
       case $host in
       *-linux-gnu)
           AC_CHECK_PROG(CXX, pgCC, pgCC)

           if test "${CXX}" = pgCC ; then 
               CC='pgcc'
               AC_DRACO_PGCC
           else
               AC_MSG_ERROR("Did not find PGI C++ compiler!")
           fi
       ;;
       *)
           AC_MSG_ERROR("PGI only available on LINUX.")
       ;;
       esac        

   elif test "${with_cxx}" = ibm ; then 
       AC_CHECK_PROG(CXX, xlC, xlC)
       AC_CHECK_PROG(CC, xlc, xlc)

       if test "${CXX}" = xlC ; then
	   AC_DRACO_IBM_VISUAL_AGE
       else
	   AC_MSG_ERROR("Did not find IBM Visual Age xlC compiler!")
       fi

   elif test "${with_cxx}" = asciwhite ; then 
       # asci white uses different executables depending upon
       # the mpi setup; so we check to see if mpi is on 
       # and set the executable appropriately 

       # mpi is on, use newmpxlC
       if test -n "${vendor_mpi}" && test "${with_mpi}" = vendor; then
	   AC_CHECK_PROG(CXX, newmpxlC, newmpxlC)
	   AC_CHECK_PROG(CC, newmpxlc, newmpxlc)

       # scalar build, use newxlC
       else
	   AC_CHECK_PROG(CXX, newxlC, newxlC)
	   AC_CHECK_PROG(CC, newxlc, newxlc)

       fi

       # check to make sure compiler is valid
       if test "${CXX}" = newxlC || test "${CXX}" = newmpxlC ; then
	   AC_DRACO_IBM_VISUAL_AGE
       else
	   AC_MSG_ERROR("Did not find ASCI White new(mp)xlC compiler!")
       fi

   else
       AC_MSG_ERROR("invalid compiler specification ${with_cxx}")

   fi

   # set the language to C++
   AC_LANG(C++)

])

dnl-------------------------------------------------------------------------dnl
dnl SGI CC COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_SGI_CC, [dnl

   AC_MSG_CHECKING("configuration of ${CXX}/${CC} compilers")

   # dirs to clean
   dirstoclean='ii_files'

   # LINKER AND LIBRARY (AR)
   LD='${CXX}'
   AR='${CXX}'
   ARLIBS='${DRACO_LIBS}'
   ARTESTLIBS='${PKG_LIBS} ${DRACO_TEST_LIBS} ${DRACO_LIBS}'

   # for CC we need to add a flag to AR to determine whether we build 
   # shared or archive libraries
   if test "${enable_shared}" = yes ; then
       ARFLAGS='-shared -o'
   else
       ARFLAGS='-ar -o'
   fi

   # COMPILATION FLAGS

   # strict asci compliance
   if test "${enable_strict_ansi:=yes}" = yes ; then
       # not really sure what the CC strict flag is, however, since
       # KCC can do our "strict" checking for us this is probably
       # not a big deal
       STRICTFLAG=""
   fi

   # optimization level
   # as opposed to KCC, -g overrides the optimization level, thus, we
   # assume that debug is the default, however, if an optimization
   # level is set we turn of debugging
   if test "${with_opt:=0}" != 0 ; then
       CXXFLAGS="${CXXFLAGS} -O${with_opt}"
       CFLAGS="${CFLAGS} -O${with_opt}" 
       enable_debug="no"
   fi

   if test "${enable_debug:=yes}" = yes ; then
       CXXFLAGS="${CXXFLAGS} -g"
       CFLAGS="${CFLAGS} -g"
   fi

   # static linking option
   if test "${enable_static_ld}" = yes ; then
       LDFLAGS="${LDFLAGS} -non_shared"
   fi

   # final compiler additions
   CXXFLAGS="${CXXFLAGS} -LANG:std -no_auto_include"
   LDFLAGS="${LDFLAGS} -LANG:std"

   AC_MSG_RESULT("SGI CC compiler flags set")

   dnl end of AC_DRACO_CC
])

dnl-------------------------------------------------------------------------dnl
dnl GNU COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_GNU_GCC, [dnl

   # finding path of gcc compiler
   AC_PATH_PROG(GCC_BIN, g++, null)

   AC_MSG_CHECKING("Setting library path of GNU compiler")
   if test "${GCC_BIN}" = null ; then
       GCC_LIB_DIR='/usr/lib'
   else
       GCC_BIN=`dirname ${GCC_BIN}`
       GCC_HOME=`dirname ${GCC_BIN}`

       # Ensure that libraries exist at this location.  If we can't
       # find libstdc++.a at this location we leave GCC_LIB_DIR set to
       # null and issue a warning.

       if test -r ${GCC_HOME}/lib/libstdc++.a; then
         GCC_LIB_DIR="${GCC_HOME}/lib"
       fi
   fi
   AC_MSG_RESULT("${GCC_LIB_DIR}")

   if test -z ${GCC_LIB_DIR}; then
       AC_MSG_WARN("Could not determine location of gcc libraries. GCC_LIB_DIR is null")
   fi

   # do compiler configuration
   AC_MSG_CHECKING("configuration of ${CXX}/${CC} compilers")

   # LINKER AND LIBRARY (AR)
   LD='${CXX}'

   # if shared then ar is gcc
   if test "${enable_shared}" = yes ; then
       AR="${CXX}"
       ARFLAGS='-shared -o'
   else
       AR='ar'
       ARFLAGS='cr'
   fi

   ARLIBS=''
   ARTESTLIBS=''

   # COMPILATION FLAGS

   # strict asci compliance
   if test "${enable_strict_ansi:=yes}" = yes ; then
       STRICTFLAG="-ansi -Wnon-virtual-dtor -Wreturn-type -pedantic"
   fi

   # optimization level
   # gcc allows -g with -O (like KCC)

   # set opt level in flags
   gcc_opt_flags="-O${with_opt:=0}"

   # set up compiler when optimized
   if test "${with_opt}" != 0; then

       # set up inlining when optimization is on
       gcc_opt_flags="-finline-functions ${gcc_opt_flags}"

       # turn off debug flag by default if not requested explicitly
       if test "${enable_debug:=no}" = yes ; then
	   gcc_opt_flags="-g ${gcc_opt_flags}"
       fi

   # set up compiler when not optimized
   else

       # default is to have debug flag on when opt=0
       if test "${enable_debug:=yes}" = yes ; then
	   gcc_opt_flags="-g ${gcc_opt_flags}"
       fi

   fi

   # add opt flags
   CXXFLAGS="${gcc_opt_flags} ${CXXFLAGS}"
   CFLAGS="${gcc_opt_flags} ${CFLAGS}"

   # RPATH FLAGS

   # add -rpath for the compiler library (G++ as LD does not do this
   # automatically) if required.
   case $host in

   # Darwin doesn't need any special flags
   powerpc-apple-darwin*)
   ;;

   # COMPAQ -> CXX
   alpha*-dec-osf*)
   ;;

   # EVERYTHING ELSE -> linux?
   *)
      if test -n "${GCC_LIB_DIR}"; then
           RPATH="${RPATH} -Xlinker -rpath ${GCC_LIB_DIR}"
      fi
   ;;
   esac

   # static linking option
   if test "${enable_static_ld}" = yes ; then
       LDFLAGS="${LDFLAGS} -Bstatic"
   fi

   AC_MSG_RESULT("GNU g++ compiler flags set")

   dnl end of AC_DRACO_GNU_GCC
])

dnl-------------------------------------------------------------------------dnl
dnl PGI COMPILER SETUP
dnl 
dnl Note that this implementation of PGI uses options that are only
dnl valid for LINUX
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_PGCC, [dnl

   # do compiler configuration
   AC_MSG_CHECKING("configuration of ${CXX}/${CC} compilers")

   # LINKER AND LIBRARY (AR)
   LD='${CXX}'

   # if shared then ar is pgCC
   if test "${enable_shared}" = yes ; then
       AR="${CXX}"
       ARFLAGS='-shared -o'

       # must use position-independent code
       CXXFLAGS="${CXXFLAGS} -fPIC"
       CFLAGS="${CFLAGS} -fPIC"
   else
       AR='ar'
       ARFLAGS='cr'
   fi

   ARLIBS=''
   ARTESTLIBS=''

   # COMPILATION FLAGS

   # strict asci compliance
   if test "${enable_strict_ansi:=yes}" = yes ; then
       STRICTFLAG="-Xa -A --no_using_std"

       # suppress long long errors in the platform-dependent options
       # section 

       # suppress missing return statement warning (we get this in
       # nearly every STL inclusion through PGICC)
       STRICTFLAG="--diag_suppress 940 ${STRICTFLAG}"
   fi

   # optimization level
   # pgCC allows -g with -O

   # set opt level in flags
   pgcc_opt_flags="-O${with_opt:=0}"

   # set up compiler when optimized
   if test "${with_opt}" != 0; then

       # set up inlining when optimization is on
       pgcc_opt_flags="${pgcc_opt_flags}"

       # turn off debug flag by default if not requested explicitly
       if test "${enable_debug:=no}" = yes ; then
	   pgcc_opt_flags="-g ${pgcc_opt_flags}"
       fi

   # set up compiler when not optimized
   else

       # default is to have debug flag on when opt=0
       if test "${enable_debug:=yes}" = yes ; then
	   pgcc_opt_flags="-g ${pgcc_opt_flags}"
       fi

   fi

   # add opt flags
   CXXFLAGS="${pgcc_opt_flags} ${CXXFLAGS}"
   CFLAGS="${pgcc_opt_flags} ${CFLAGS}"
   
   # add ieee flag
   CXXFLAGS="${CXXFLAGS} -Kieee"
   CFLAGS="${CFLAGS} -Kieee"

   # instantiate only functions that are used in the compilation
   CXXFLAGS="${CXXFLAGS} -t=used --no_implicit_include"

   # set unnormalized values to zero
   CXXFLAGS="${CXXFLAGS} -Mdaz"
   CFLAGS="${CFLAGS} -Mdaz"

   AC_MSG_RESULT("PGI pgCC compiler flags set")

   dnl end of AC_DRACO_PGCC
])

dnl-------------------------------------------------------------------------dnl
dnl COMPAQ CXX COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_COMPAQ_CXX, [dnl

   dnl 6-FEB-02 NEED TO ADD MODS !!!!!!

   AC_MSG_CHECKING("configuration of ${CXX}/${CC} compilers")

   # CXX SPECIFIC FLAGS
   dirstoclean='cxx_repository'

   # LINKER AND LIBRARY (AR)
   LD='${CXX}'

   # if shared then ar is cxx
   if test "${enable_shared}" = yes ; then
       AR='${CXX}'
       ARFLAGS="-shared -nocxxstd"
       ARFLAGS="${ARFLAGS} -o"
   else
       AR='ar'
       ARFLAGS='cr'
   fi

   # the contents of the cxx_repository do not seem to need adding 
   # when building shared libraries; you do have to add them for
   # archives 
   if test "${enable_shared}" != yes ; then
       ARLIBS='$(wildcard cxx_repository/*)'
       ARTESTLIBS='$(wildcard cxx_repository/*)'
   fi

   # COMPILATION FLAGS

   # strict asci compliance
   if test "${enable_strict_ansi:=yes}" = yes ; then
       STRICTFLAG="-std strict_ansi"
       CXX="${CXX} -model ansi"
   fi

   # make sure we always use the standard IO stream
   CPPFLAGS="${CPPFLAGS} -D__USE_STD_IOSTREAM" 

   # optimization level

   # if optimization is on turn off debug flag unless asked for
   if test "${with_opt:=0}" != 0 ; then

       # if debug is on then use -g1,2,3
       if test "${enable_debug:=no}" = yes ; then
	   cxx_opt_flag="-g${with_opt}"
       else
	   cxx_opt_flag="-O${with_opt}"
       fi

   # turn off optimizations
   else
   
       # we want -g unless not asked for
       if test "${enable_debug:=yes}" = yes ; then
	   cxx_opt_flag="-g -O0"
       else
	   cxx_opt_flag="-O0"
       fi

   fi

   # set up cxx flags
   CXXFLAGS="${CXXFLAGS} ${cxx_opt_flag}"
   CFLAGS="${CFLAGS} ${cxx_opt_flag}"

   # add ieee flag
   CXXFLAGS="${CXXFLAGS} -ieee"
   CFLAGS="${CFLAGS} -ieee"

   # turn off implicit inclusion
   CXXFLAGS="${CXXFLAGS} -noimplicit_include"

   # use the -pt template option for the compiler:
   # -pt Automatically instantiate templates into the repository with
   #  external linkage. Manually instantiated templates are placed in
   #  the output object with external linkage. This option is the default.
   CXXFLAGS="${CXXFLAGS} -pt"

   # static linking option
   if test "${enable_static_ld}" = yes ; then
       LDFLAGS="${LDFLAGS} -non_shared"
   fi

   # add thread safe linkage
   LDFLAGS="${LDFLAGS}" # -pthread"

   AC_MSG_RESULT("CXX Compaq compiler flags set")
   
   dnl end of AC_DRACO_COMPAQ_CXX
])

dnl-------------------------------------------------------------------------dnl
dnl Intel icpc COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_INTEL_ICPC, [dnl

   AC_MSG_CHECKING("configuration of ${CXX}/${CC} compilers")

   # icpc SPECIFIC FLAGS

   # LINKER AND LIBRARY
   LD='${CXX}'

   # if shared then ar is icpc
   if test "${enable_shared}" = yes ; then
       AR="${CXX}"
       ARFLAGS='-shared -o'
   else
       AR='ar'
       ARFLAGS='cr'
   fi

   ARLIBS=''
   ARTESTLIBS=''

   # COMPILATION FLAGS

   # strict asci compliance
   if test "${enable_strict_ansi:=yes}" = yes ; then
       STRICTFLAG="-ansi"
   fi

   # set up compiler when optimized (enable inline keyword but not
   # compiler-choice inlining)
   if test "${with_opt:=0}" != 0 ; then

       # turn off debug by default
       if test "${enable_debug:=no}" = yes ; then
	   icpc_opt_flags="-g -O${with_opt} -Ob1 -ip"
       else
	   icpc_opt_flags="-O${with_opt} -Ob1"
       fi

   #set up compiler when not optimized (turn off inlining with -Ob0)
   else

       # turn on debug by default
       if test "${enable_debug:=yes}" = yes ; then
	   icpc_opt_flags="-g -O0 -Ob0"
       else
	   icpc_opt_flags="-O0 -Ob0"
       fi

   fi
   
   # set the cxx and c flags
   CXXFLAGS="${CXXFLAGS} ${icpc_opt_flags}"
   CFLAGS="${CFLAGS} ${icpc_opt_flags}"

   # static linking option
   if test "${enable_static_ld}" = yes ; then
       LDFLAGS="${LDFLAGS} -static"
   fi

   AC_MSG_RESULT("icpc compiler flags set")
   
   dnl end of AC_DRACO_INTEL_ICPC
])

dnl-------------------------------------------------------------------------dnl
dnl IBM VISUAL AGE COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_IBM_VISUAL_AGE, [dnl

   AC_MSG_CHECKING("configuration of ${CXX}/${CC} compilers")

   # xlC SPECIFIC FLAGS

   # LINKER AND LIBRARY
   LD='${CXX}'

   # if shared then ar is xlC
   if test "${enable_shared}" = yes ; then
       AR="${CXX}"
       ARFLAGS='-brtl -Wl,-bh:5 -G -o'

       # when AR=newmpxlC we need to add /lib/crt0.o to 
       # avoid p_argcx and p_argvx link error when building libs
       if test "${AR}" = newmpxlC ; then
	   ARLIBS='/lib/crt0.o'
	   ARTESTLIBS='/lib/crt0.o'
       fi

       ARLIBS="${ARLIBS} \${DRACO_LIBS} \${VENDOR_LIBS}"
       ARTESTLIBS="${ARTESTLIBS} \${PKG_LIBS} \${DRACO_TEST_LIBS}"
       ARTESTLIBS="${ARTESTLIBS} \${DRACO_LIBS}\${VENDOR_TEST_LIBS}"
       ARTESTLIBS="${ARTESTLIBS} \${VENDOR_LIBS}"
   else
       AR='ar'
       ARFLAGS='cr'

       ARLIBS=''
       ARTESTLIBS=''
   fi

   # COMPILATION FLAGS

   # strict asci compliance
   if test "${enable_strict_ansi:=yes}" = yes ; then
       STRICTFLAG="-qlanglvl=strict98"
   fi

   # the qinline option controls inlining, when -g is on no inlining
   # is done, with -O# inlining is on by default

   # set up compiler when optimized 
   if test "${with_opt:=0}" != 0; then

       # optflags
       xlC_opt_flags="-qarch=auto -qtune=auto -qcache=auto"

       # optimization level    
       if test "${with_opt}" = 1; then
	   # if asking for 1 just use opt in ibm   
	   xlC_opt_flags="${xlC_opt_flags} -qopt"
       else
	   # otherwise use number

	   # turn of aggressive semantic optimizations on all levels
	   # -O2 and above
	   xlC_opt_flags="${xlC_opt_flags} -qopt=${with_opt} -qstrict"
       fi

       # turn off debug by default
       if test "${enable_debug:=no}" = yes ; then
	   xlC_opt_flags="-g ${xlC_opt_flags}"
       fi

   #set up compiler when not optimized 
   else

       # optflags
       xlC_opt_flags="-qnoopt"

       # turn on debug by default
       if test "${enable_debug:=yes}" = yes ; then
	   xlC_opt_flags="-g ${xlC_opt_flags}"
       fi

   fi
   
   # set the CXX and CC flags

   # set the optimizations
   CXXFLAGS="${CXXFLAGS} ${xlC_opt_flags}"
   CFLAGS="${CFLAGS} ${xlC_opt_flags}"

   # set template stuff
   CXXFLAGS="${CXXFLAGS} -w -qnotempinc"

   # static linking option
   if test "${enable_static_ld:=no}" = yes ; then
       LDFLAGS="${LDFLAGS} -bstatic"

   # if we are building shared libraries we need to add
   # run-time-linking
   else
       LDFLAGS="${LDFLAGS} -brtl -Wl,-bh:5"

   fi

   AC_MSG_RESULT("${CXX} compiler flags set")
   
   dnl end of AC_DRACO_IBM_VISUAL_AGE
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_compiler.m4
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl ac_platforms.m4
dnl
dnl Defines platform-specfic environments, including default vendor
dnl settings for the CCS-4/ASC computer platforms.
dnl
dnl Thomas M. Evans
dnl 2003/04/30 20:29:39
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_PLATFORM_ENVIRONMENT
dnl
dnl Configure draco build system platfrom-specfic variables
dnl This function is called within AC_DRACO_ENV
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_PLATFORM_ENVIRONMENT], [dnl

   # we must know the host
   AC_REQUIRE([AC_CANONICAL_HOST])

   # dependency rules
   DEPENDENCY_RULES='Makefile.dep.general'

   # systems setup
   case $host in

   # ***********
   # LINUX SETUP
   # ***********
   *-linux-gnu)
       AC_DBS_LINUX_ENVIRONMENT
   ;;

   # ***********
   # CYGWIN SETUP
   # ***********
   i686-pc-cygwin)
       AC_DBS_CYGWIN_ENVIRONMENT
   ;;

   # *********
   # SGI SETUP
   # *********
   mips-sgi-irix6.*)
       AC_DBS_IRIX_ENVIRONMENT
   ;;

   # ******************
   # TRU64 COMPAQ SETUP
   # ******************
   alpha*-dec-osf*)
       AC_DBS_OSF_ENVIRONMENT
   ;;

   # *************
   # IBM AIX SETUP
   # *************
   *ibm-aix*)
       AC_DBS_IBM_ENVIRONMENT
   ;;

   # *****************
   # SUN/SOLARIS SETUP
   # *****************
   sparc-sun-solaris2.*)
       AC_DBS_SUN_ENVIRONMENT
   ;;

   # *********************
   # MAC OS X/DARWIN SETUP
   # *********************
   powerpc-apple-darwin*)
       AC_DBS_DARWIN_ENVIRONMENT
   ;;      

   # *******
   # NOTHING
   # *******
   *)
       AC_MSG_ERROR("Cannot figure out the platform or host!")
   ;;
   esac
])

dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_LAHEY_ENVIRONMENT
dnl
dnl Some vendor setups require that the Lahey lib dir and compiler
dnl libraries be provided on the link line.  This m4 function adds the
dnl necessary libraries to LIBS.
dnl-------------------------------------------------------------------------dnl
AC_DEFUN([AC_DBS_LAHEY_ENVIRONMENT], [dnl

   if test "${CXX}" != g++; then
       AC_MSG_ERROR("LAHEY must be configured with g++ on LINUX.")
   fi

   AC_MSG_CHECKING("for extra lf95 library requirements.")
   if test -n "${vendor_eospac}"    ||
      test -n "${vendor_scalapack}" ||
      test -n "${vendor_trilinos}"; then
         f90_lib_loc=`which lf95 | sed -e 's/bin\/lf95/lib/'`
	 extra_f90_libs="-L${f90_lib_loc} -lfj9i6 -lfj9e6 -lfj9f6 -lfst -lfccx86_6a"
         LIBS="${LIBS} ${extra_f90_libs}"
         extra_f90_rpaths="-Xlinker -rpath ${f90_lib_loc}"
         RPATH="${RPATH} ${extra_f90_rpaths}"
         AC_MSG_RESULT("${extra_f90_libs}")
   else
         AC_MSG_RESULT("none.")
   fi

   dnl Optimize flag   
   AC_MSG_CHECKING("for F90FLAGS")
   if test "${with_opt:=0}" != 0 ; then
      F90FLAGS="${F90FLAGS} -O${with_opt}"
   else 
      F90FLAGS="${F90FLAGS} -g"
   fi

   dnl C preprocessor flag
   F90FLAGS="${F90FLAGS} -Cpp"
   AC_MSG_RESULT(${F90FLAGS})

   dnl scalar or mpi ?
   AC_MSG_CHECKING("for F90MPIFLAGS")
   if test ${with_mpi:=no} = "no"; then
      F90MPIFLAGS="${F90FLAGS} -DC4_SCALAR"
   else
      if test "${with_mpi}" = mpich; then
         F90MPIFLAGS="-lfmpich"
      else dnl LAMPI support
         F90MPIFLAGS="-lmpi"
      fi
   fi
   AC_MSG_RESULT(${F90MPIFLAGS})
   
   dnl Add C++ options to F90 link line
   AC_MSG_CHECKING("for F90CXXFLAGS")
   CXXLIBDIR=${GCC_LIB_DIR}
   F90CXXFLAGS="-L${CXXLIBDIR} -lstdc++"
   AC_MSG_RESULT(${F90CXXFLAGS})

   AC_MSG_CHECKING("for F90VENDOR_LIBS")
   F90VENDOR_LIBS="$F90VENDOR_LIBS ${F90MPIFLAGS} ${F90CXXFLAGS}"
   AC_MSG_RESULT("${F90VENDOR_LIBS}")
])

dnl ------------------------------------------------------------------------ dnl
dnl AC_DBS_PGF90_ENVIRONMENT
dnl
dnl Some vendor setups require that the Portland Group F90 lib dir and
dnl compiler libraries be provided on the link line.  This m4 function
dnl adds the necessary libraries to LIBS.
dnl ------------------------------------------------------------------------ dnl
AC_DEFUN([AC_DBS_PGF90_ENVIRONMENT], [dnl

   # set the proper RPATH command depending on the C++ compiler
   case ${CXX} in 
       g++ | icpc)
           rpath='-Xlinker -rpath '
           ;;
       pgCC)
           rpath='-R'
           ;;
       *)
           AC_MSG_ERROR("Improper compiler set in LINUX.")
   esac

   AC_MSG_CHECKING("for extra pgf90 library requirements.")
   if test -n "${vendor_eospac}"    ||
      (test -n "${vendor_lapack}" && test "${with_lapack}" = "atlas") ||
      test -n "${vendor_scalapack}" ||
      test -n "${vendor_trilinos}"; then
         f90_lib_loc=`which pgf90 | sed -e 's/bin\/pgf90/lib/'`
         if test -r ${f90_lib_loc}/libpgc.a; then
	    extra_f90_libs="-L${f90_lib_loc} -lpgf90 -lpgf902 -lpgc -lpgftnrtl"
            extra_f90_libs="${extra_f90_libs} -lpgf90_rpm1 -lpghpf2"
            extra_f90_rpaths="$rpath${f90_lib_loc}"
         else
	    extra_f90_libs="-L${f90_lib_loc} -lpgf90 -lpgf902 -lpgftnrtl"
            extra_f90_libs="${extra_f90_libs} -lpgf90_rpm1 -lpghpf2"
            f90_lib_loc2=`which pgf90 | sed -e 's/bin\/pgf90/lib-linux86-g232/'`
            if test -r ${f90_lib_loc2}/libpgc.a; then
               extra_f90_libs="${extra_f90_libs} -L${f90_lib_loc2} -lpgc"
               extra_f90_rpaths="-Xlinker -rpath ${f90_lib_loc}"
               extra_f90_rpaths="${extra_f90_rpaths} $rpath${f90_lib_locs}"
            fi
         fi
         LIBS="${LIBS} ${extra_f90_libs}"
         RPATH="${RPATH} ${extra_f90_rpaths}"
         AC_MSG_RESULT("${extra_f90_libs}")
   else
         AC_MSG_RESULT("none.")
   fi

   dnl Optimize flag   
   AC_MSG_CHECKING("for F90FLAGS")
   if test "${with_opt:=0}" != 0 ; then
      if test ${with_opt} -gt 2; then
         F90FLAGS="${F90FLAGS} -O2"
      else
         F90FLAGS="${F90FLAGS} -O${with_opt}"
      fi
   else 
      F90FLAGS="${F90FLAGS} -g"
   fi

   dnl C preprocessor flag
   F90FLAGS="${F90FLAGS} -Mpreprocess"
   AC_MSG_RESULT(${F90FLAGS})

   dnl scalar or mpi ?
   AC_MSG_CHECKING("for F90MPIFLAGS")
   if test ${with_mpi:=no} = "no"; then
      F90FLAGS="${F90FLAGS} -DC4_SCALAR"
   else
      if test "${with_mpi}" = mpich; then
         F90MPIFLAGS="-lfmpich"
      else dnl LAMPI support
         F90MPIFLAGS="-lmpi"
      fi
   fi
   AC_MSG_RESULT(${F90MPIFLAGS})

   dnl Add C++ options to F90 link line
   AC_MSG_CHECKING("for F90CXXFLAGS")
   CXXLIBDIR=${GCC_LIB_DIR}
   F90CXXFLAGS="-L${CXXLIBDIR} -lstdc++"
   AC_MSG_RESULT(${F90CXXFLAGS})

   AC_MSG_CHECKING("for F90VENDOR_LIBS")
   F90VENDOR_LIBS="$F90VENDOR_LIBS ${F90MPIFLAGS} ${F90CXXFLAGS}"
   AC_MSG_RESULT("${F90VENDOR_LIBS}")
])

dnl ------------------------------------------------------------------------ dnl
dnl AC_DBS_COMPAQ_F90_ENVIRONMENT
dnl
dnl Some vendor setups require that the Portland Group F90 lib dir and
dnl compiler libraries be provided on the link line.  This m4 function
dnl adds the necessary libraries to LIBS.
dnl ------------------------------------------------------------------------ dnl
AC_DEFUN([AC_DBS_COMPAQ_F90_ENVIRONMENT], [dnl

   f90_lib_loc=`which f90 | sed -e 's/bin\/f90/lib/'`
   cxx_lib_loc=`which cxx | sed -e 's/bin\/cxx/lib/'`

   AC_MSG_CHECKING("for extra f90 library requirements.")
   if test -n "${vendor_eospac}"    ||
      test -n "${vendor_gandolf}"   || 
      test -n "${vendor_pcg}"       || 
      test -n "${vendor_udm}"       ||
      test -n "${vendor_blacs}"; then

      extra_f90_libs="-L${f90_lib_loc} -lfor"
      extra_f90_rpaths="-rpath ${f90_lib_loc}"

      LIBS="${LIBS} ${extra_f90_libs}"
      RPATH="${RPATH} ${extra_f90_rpaths}"
      AC_MSG_RESULT("${extra_f90_libs}")

   else
         AC_MSG_RESULT("none.")
   fi

   dnl Optimize flag   
   AC_MSG_CHECKING("for F90FLAGS")
   if test "${with_opt:=0}" != 0 ; then
      if test ${with_opt} -gt 2; then
         F90FLAGS="${F90FLAGS} -O2"
      else
         F90FLAGS="${F90FLAGS} -O${with_opt}"
      fi
   else 
      F90FLAGS="${F90FLAGS} -g"
   fi

   dnl C preprocessor flag
   F90FLAGS="${F90FLAGS} -cpp"
   AC_MSG_RESULT(${F90FLAGS})

   dnl scalar or mpi ?
   AC_MSG_CHECKING("for F90MPIFLAGS")
   if test ${with_mpi:=no} = "no"; then
      F90FLAGS="${F90FLAGS} -DC4_SCALAR"
   else
      F90MPIFLAGS="${F90FLAGS} -lfmpi"
   fi
   AC_MSG_RESULT(${F90MPIFLAGS})

   if test -n "${vendor_pcg}"  || 
      test "${with_udm}" = mpi || 
      test -n "${vendor_blacs}" ; then
      LIBS="${LIBS} ${F90MPIFLAGS}"
   fi

   dnl Add C++ options to F90 link line
   AC_MSG_CHECKING("for F90CXXFLAGS")
   CXXLIBDIR=${cxx_lib_loc}
   F90CXXFLAGS="-L${CXXLIBDIR} -lcxxstdma -lcxxma"
   AC_MSG_RESULT(${F90CXXFLAGS})

   AC_MSG_CHECKING("for F90VENDOR_LIBS")
   F90VENDOR_LIBS="$F90VENDOR_LIBS ${F90MPIFLAGS} ${F90CXXFLAGS}"
   AC_MSG_RESULT("${F90VENDOR_LIBS}")
])

dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_LINUX_ENVIRONMENT
dnl
dnl Configure draco build system Linux-specific variables
dnl This function is called within AC_DBS_PLATFORM_ENVIRONMENT
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_LINUX_ENVIRONMENT], [dnl

       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       AC_DBS_SETUP_POSIX

       #   
       # LONG LONG on Linux
       #
       
       # always allow long long in strict ansi mode (if possible)
       
       if test -n "${STRICTFLAG}"; then

           case $CXX in

           # GNU g++
           g++) 
               AC_MSG_NOTICE([g++ -ansi option set to allow long long type!])
               STRICTFLAG="$STRICTFLAG -Wno-long-long"
           ;;

           # PGI
           pgCC)
               AC_MSG_NOTICE([pgCC supressing error 450 for long long!])
               STRICTFLAG="--diag_suppress 450 ${STRICTFLAG}"  
           ;;

           # catchall
           *) 
               # do nothing
           ;;

           esac

       fi

       # 
       # end of LONG LONG setup
       #

       #
       # Setup communications packages
       #
       AC_DBS_SETUP_COMM(mpich)

       # 
       # setup lapack 
       #
       
       # we assume that the vendor option on linux is the install of
       # redhat rpms in /usr/lib; we don't worry about atlas because
       # that has already been defined

       if test "${with_lapack}" = vendor ; then
	   lapack_libs='-llapack -lblas'
       fi 

       # 
       # end of lapack setup
       # 

       # setup F90 libs, rpath, etc.
       AC_CHECK_PROGS(F90, pgf90 lf95)
       case ${F90} in
       lf95)
          AC_DBS_LAHEY_ENVIRONMENT
          ;;
       pgf90)
          AC_DBS_PGF90_ENVIRONMENT
          ;;
       esac

       #
       # add libg2c to LIBS if lapack, gandolf, or pcg is used
       #
       AC_MSG_CHECKING("libg2c requirements")
       if test -n "${vendor_lapack}" || 
          test -n "${vendor_pcg}"    ||
	  test -n "${vendor_gandolf}"; then

	   # Add g2c for various compilers
           case "${CXX}" in

               g++)
                   LIBS="${LIBS} -lg2c"
                   AC_MSG_RESULT("-lg2c added to LIBS")
                   ;;

               icpc | pgCC)
                   AC_PATH_PROG(GCC_BIN, g++, null)
                   GCC_BIN=`dirname ${GCC_BIN}`
                   GCC_HOME=`dirname ${GCC_BIN}`
                   GCC_LIB_DIR="${GCC_HOME}/lib"
                   LIBS="${LIBS} -L${GCC_LIB_DIR} -lg2c"
                   AC_MSG_RESULT("-lg2c added to LIBS")
                   ;;

               *)
                   AC_MSG_RESULT("not needed")
                   ;;

           esac

       else
	   AC_MSG_RESULT("not needed")
       fi

       #
       # add librt to LIBS if udm is used
       #
       AC_MSG_CHECKING("librt requirements")
       if test -n "${vendor_udm}"; then

	   # Add rt for g++
	   if test "${CXX}" = g++ ; then
	       LIBS="${LIBS} -lrt"
	       AC_MSG_RESULT("-lrt added to LIBS")
           else
               AC_MSG_RESULT("not needed")
	   fi

       else
           AC_MSG_RESULT("not needed")
       fi

       #
       # If dlopen is specified, 1) add libdl to LIBS; 
       # 2) add -fPIC to compile flags.
       #
       AC_MSG_CHECKING("libdl requirements")
       if test -n "${vendor_dlopen}" ; then
           if test "${enable_dlopen}" = yes ; then
               LIBS="${LIBS} -ldl"

               # if we are using g++ add fPIC (pgCC already has fPIC
               # when building shared libraries
               if test "${CXX}" = g++; then
                   CFLAGS="${CFLAGS} -fPIC"
                   CXXFLAGS="${CXXFLAGS} -fPIC"
                   AC_MSG_RESULT("-ldl added to LIBS -fPIC added to compile flags")
               else
                   AC_MSG_RESULT("-ldl added to LIBS")
               fi

           else  
               AC_MSG_RESULT("not needed")
           fi

       else
           AC_MSG_RESULT("not needed")
       fi

       #
       # Set up fpe_trap for this platform if gcc is on.
       #
       if test "${CXX}" = g++; then
           AC_DEFINE(FPETRAP_LINUX_X86)
       fi

       #
       # finalize vendors
       #
       AC_VENDOR_FINALIZE

       # handle rpaths
       case ${CXX} in
           pgCC)
               AC_DBS_SETUP_RPATH(-R, nospace)
               ;;
           g++ | icpc)
               AC_DBS_SETUP_RPATH('-Xlinker -rpath', space)
               ;;
           *)
               AC_MSG_ERROR("Unrecognized compiler on LINUX")
               ;;
       esac

       # add the intel math library for better performance when
       # compiling with intel
       if test "${CXX}" = icpc; then
	   LIBS="$LIBS -limf"
       fi
]) dnl linux


dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_CYWGIN_ENVIRONMENT
dnl
dnl Configure draco build system Cygwin-specific variables
dnl This function is called within AC_DBS_PLATFORM_ENVIRONMENT
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_CYGWIN_ENVIRONMENT], [dnl

       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       AC_DBS_SETUP_POSIX

       AC_DBS_SETUP_COMM(mpich)

       dnl 
       dnl setup lapack 
       dnl
       
       dnl we assume that the vendor option on linux is the install of
       dnl redhat rpms in /usr/lib; we don't worry about atlas because
       dnl that has already been defined

       if test "${with_lapack}" = vendor ; then
	   lapack_libs='-llapack -lblas'
       fi 

       dnl 
       dnl end of lapack setup
       dnl 

       # setup lf95 libs
       AC_DBS_LAHEY_ENVIRONMENT

       dnl
       dnl add libg2c to LIBS if lapack, gandolf, or pcg is used
       dnl
       AC_MSG_CHECKING("libg2c requirements")
       if test -n "${vendor_lapack}" || test -n "${vendor_pcg}" ||
	  test -n "${vendor_gandolf}"; then
	   
	   dnl Add g2c for various compilers
	   if test "${CXX}" = g++ ; then
	       LIBS="${LIBS} -lg2c"
	       AC_MSG_RESULT("-lg2c added to LIBS")
	   fi

       else
	   AC_MSG_RESULT("not needed")
       fi

       dnl
       dnl finalize vendors
       dnl
       AC_VENDOR_FINALIZE

       AC_DBS_SETUP_RPATH('-Xlinker -rpath', space)

]) dnl cygwin

dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_OSF_ENVIRONMENT
dnl
dnl Configure draco build system OSF-specific variables
dnl This function is called within AC_DBS_PLATFORM_ENVIRONMENT
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_OSF_ENVIRONMENT], [dnl

       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       AC_DBS_SETUP_POSIX

       #
       # setup communication packages
       #

       # setup vendor mpi
       if test "${with_mpi}" = vendor ; then

	   # define mpi libraries
	   # note: mpi and mpio are separate libraries on compaqs
	   mpi_libs='-lmpi -lmpio'
       
       # setup mpich
       elif test "${with_mpi}" = mpich ; then

	   # define mpi libraries
	   mpi_libs='-lmpich'
   
       fi

       # add COMPAQ ALASKA Specfic options
       if test "${with_mpi}" = vendor ; then
	   # define version check
	   AC_DEFINE(MPI_VERSION_CHECK)
       fi

       #
       # end of communication packages
       #

       #
       # setup lapack
       #

       AC_MSG_CHECKING("for lapack libraries")
       if test "${with_lapack}" = vendor ; then
	   lapack_libs='-lcxmlp -lcxml'
           AC_MSG_RESULT("${lapack_libs}")
       else
           AC_MSG_RESULT("none.")
       fi

       #
       # end of lapack setup
       #

       #
       # udm requires long long warnings to be disabled
       #

       if test -n "${vendor_udm}" ; then
	   STRICTFLAG="${STRICTFLAG} -msg_disable nostdlonglong"
	   STRICTFLAG="${STRICTFLAG} -msg_disable nostdlonglong"
       fi

       #
       # end of udm setup
       #

       #
       # FORTRAN configuration for Compaq f90
       # setup F90, libs, rpath, etc.
       #
       AC_CHECK_PROGS(F90, f90)
       case ${F90} in
       f90)
          AC_DBS_COMPAQ_F90_ENVIRONMENT
          ;;
       esac

       #
       # libudm/librmscall setup
       #

       AC_MSG_CHECKING("librmscall requirements")
       if test -n "${vendor_udm}"; then
          LIBS="${LIBS} -lrmscall"
          AC_MSG_RESULT("-lrmscall added to LIBS")
       else
	   AC_MSG_RESULT("not needed")
       fi

       #
       # end of libudm setup
       #

       #
       # Set up fpe_trap for this platform.
       #
       AC_DEFINE(FPETRAP_OSF_ALPHA)

       #
       # finalize vendors
       #
       AC_VENDOR_FINALIZE

       AC_DBS_SETUP_RPATH(-rpath, colon)

]) dnl osf

dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_IBM_ENVIRONMENT
dnl
dnl Configure draco build system IBM-specific variables
dnl This function is called within AC_DBS_PLATFORM_ENVIRONMENT
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_IBM_ENVIRONMENT], [dnl

       # dependency rules for IBM visual age compiler are complex
       if test "${with_cxx}" = asciwhite || test "${with_cxx}" = ibm; then
	   DEPENDENCY_RULES='Makefile.dep.xlC'
       fi
   
       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       AC_DBS_SETUP_POSIX

       # set up 32 or 64 bit compiling on IBM
       if test "${enable_32_bit:=no}" = yes ; then
	   
	   # switch on gcc or xlC compiler
	   if test "${with_cxx}" = gcc; then
	       CXXFLAGS="${CXXFLAGS} -maix32"
	       CFLAGS="${CFLAGS} -maix32"
	   elif test "${with_cxx}" = asciwhite || 
                test "${with_cxx}" = ibm; then
	       CXXFLAGS="${CXXFLAGS} -q32"
	       CFLAGS="${CFLAGS} -q32"
	   fi

       elif test "${enable_64_bit:=no}" = yes ; then
	   
	   # switch on gcc or xlC compiler
	   if test "${with_cxx}" = gcc; then
	       CXXFLAGS="${CXXFLAGS} -maix64"
	       CFLAGS="${CFLAGS} -maix64"
	   elif test "${with_cxx}" = asciwhite || 
                test "${with_cxx}" = ibm; then
	       CXXFLAGS="${CXXFLAGS} -q64"
	       CFLAGS="${CFLAGS} -q64"
	   fi

       fi

       # set up the heap size
       if test "${with_cxx}" = asciwhite ; then
	   LDFLAGS="${LDFLAGS} -bmaxdata:0x80000000"
       fi

       # 
       # GCC on AIX FLAGS
       #
       if test "${with_cxx}" = gcc; then

	   # add the appropriate runtime linking for shared compiling
	   if test "${enable_shared}" = yes; then
	       ARFLAGS="-Xlinker -brtl -Xlinker -bh:5 ${ARFLAGS}"
	       ARLIBS='${DRACO_LIBS} ${VENDOR_LIBS}'
	       ARTESTLIBS='${PKG_LIBS} ${DRACO_TEST_LIBS} ${DRACO_LIBS}'
	       ARTESTLIBS="${ARTESTLIBS} \${VENDOR_TEST_LIBS} \${VENDOR_LIBS}" 
	   fi

	   # we always allow shared object linking
	   if test "${enable_static_ld}" != yes; then
	       LDFLAGS="${LDFLAGS} -Xlinker -brtl -Xlinker -bh:5"
	   fi

	   # turn off the rpath
	   RPATH=''
       fi

       #
       # setup communication packages
       #
       if test -n "${vendor_mpi}"; then

	   # setup vendor mpi
	   if test "${with_mpi}" = vendor ; then

	       # on asciwhite the newmpxlC compiler script takes care
	       # of loading the mpi libraries; since it will fail
	       # if libraries are loaded and newmpxlC is used; throw
	       # an error if it occurs
	       if test "${with_cxx}" = asciwhite; then

		   if test -n "${MPI_INC}" || test -n "${MPI_LIB}"; then
		       AC_MSG_ERROR("Cannot set mpi paths with newmpxlC.")
		   fi

		   mpi_libs=''

	       fi

	       # set up libraries if we are on ibm
	       if test "${with_cxx}" = ibm; then

		   # set up mpi library
		   mpi_libs='-lmpi'

	       fi

	       # now turn on long long support if we are using the 
	       # visual age compiler
	       if test "${with_cxx}" = ibm || 
	          test "${with_cxx}" = asciwhite ; then

		   if test "${enable_strict_ansi}"; then
		       AC_MSG_WARN("xlC set to allow long long")
		       STRICTFLAG="-qlanglvl=extended"
		       CFLAGS="${CFLAGS} -qlonglong"
		       CXXFLAGS="${CXXFLAGS} -qlonglong"
		   fi

	       fi   
       
	   # setup mpich
	   elif test "${with_mpi}" = mpich ; then

	       # set up mpi libs
	       mpi_libs='-lmpich'
   
	   fi

       fi
       #
       # end of communication packages
       #

       #
       # OTHER VENDORS
       #

       # we don't have the other vendors setup explicitly 

       #
       # finalize vendors
       #
       AC_VENDOR_FINALIZE

       # RPATH is derived from -L, we don't need an explicit setup.

       # do shared specific stuff
       if test "${enable_shared}" = yes ; then
	   # turn off ranlib
	   RANLIB=':'
       fi

])


dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_SUN_ENVIRONMENT
dnl
dnl Configure draco build system Sun-specific variables
dnl This function is called within AC_DBS_PLATFORM_ENVIRONMENT
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_SUN_ENVIRONMENT], [dnl

       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       AC_DBS_SETUP_POSIX

       #
       # setup communication packages
       #
   
       # setup for mpi support
       # we only support mpich on sgis       
       if test "${with_mpi}" = vendor ; then

	   AC_MSG_ERROR("We do not support vendor mpi on the SUN yet!")

       elif test "${with_mpi}" = mpich ; then
	   
	   # define sun-required libraries for mpich, v 1.0 (this
	   # needs to be updated for version 1.2)
	   mpi_libs='-lpmpi -lmpi -lsocket -lnsl'

       fi

       #
       # end of communication package setup
       #

       #
       # setup lapack
       #

       if test "${with_lapack}" = vendor ; then
	   lapack_libs='-llapack -lblas -lF77 -lM77 -lsunmath'
       fi

       #
       # end of lapack setup
       #

       #
       # finalize vendors
       #
       AC_VENDOR_FINALIZE

       AC_DBS_SETUP_RPATH(-R, space)
]) dnl sun

dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_IRIX_ENVIRONMENT
dnl
dnl Configure draco build system IRIX-specific variables
dnl This function is called within AC_DBS_PLATFORM_ENVIRONMENT
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_IRIX_ENVIRONMENT], [dnl

       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       AC_DBS_SETUP_POSIX

       # RANLIB TAG ON SGI
       RANLIB=':'

       # BIT COMPILER FLAGS ON SGI
       if test "${enable_32_bit:=no}" = yes ; then
	   if test "${with_cxx}" = gcc ; then
	       CXXFLAGS="-mabi=n32 ${CXXFLAGS}"
	       CFLAGS="-mabi=n32 ${CFLAGS}"
	       if test "${enable_shared}" = yes ; then
		   ARFLAGS="-mabi=n32 ${ARFLAGS}"
	       fi
	       LDFLAGS="-mabi=n32 ${LDFLAGS}"
	   else
	       CXXFLAGS="-n32 ${CXXFLAGS}"
	       CFLAGS="-n32 ${CFLAGS}"
	       ARFLAGS="-n32 ${ARFLAGS}"
	       LDFLAGS="-n32 ${LDFLAGS}"
	   fi
       else 
	   if test "${with_cxx}" = gcc ; then
	       CXXFLAGS="-mabi=64 ${CXXFLAGS}"
	       CFLAGS="-mabi=64 ${CFLAGS}"
	       if test "${enable_shared}" = yes ; then
		   ARFLAGS="-mabi=64 ${ARFLAGS}"
	       fi
	       LDFLAGS="-mabi=64 ${LDFLAGS}"
	   else
	       CXXFLAGS="-64 ${CXXFLAGS}"
	       CFLAGS="-64 ${CFLAGS}"
	       ARFLAGS="-64 ${ARFLAGS}"
	       LDFLAGS="-64 ${LDFLAGS}"
	   fi
       fi

       # MIPS INSTRUCTIONS ON SGI
       # this is different depending upon the compiler
       if test "${with_cxx}" = kcc ; then
	   CXXFLAGS="-mips${with_mips:=4} --backend -r10000 ${CXXFLAGS}"
	   CFLAGS="-mips${with_mips:=4} -r10000 ${CFLAGS}"
	   ARFLAGS="-mips${with_mips:=4} ${ARFLAGS}"
	   LDFLAGS="-mips${with_mips:=4} ${LDFLAGS}"
       elif test "${with_cxx}" = sgi ; then
	   CXXFLAGS="-mips${with_mips:=4} -r10000 ${CXXFLAGS}"
	   CFLAGS="-mips${with_mips:=4} -r10000 ${CFLAGS}"
	   ARFLAGS="-mips${with_mips:=4} ${ARFLAGS}"
	   LDFLAGS="-mips${with_mips:=4} ${LDFLAGS}"
       elif test "${with_cxx}" = gcc ; then
	   CXXFLAGS="-mips${with_mips:=4} ${CXXFLAGS}"
	   CFLAGS="-mips${with_mips:=4} ${CFLAGS}"
	   if test "${enable_shared}" = yes ; then
	       ARFLAGS="-mips${with_mips:=4} ${ARFLAGS}"
	   fi
	   LDFLAGS="-mips${with_mips:=4} ${LDFLAGS}"
       fi

       #
       # setup communication packages
       #
   
       # setup for mpi support
       # we only support vendor mpi on sgis       
       if test "${with_mpi}" = vendor ; then
	   
	   # mpi library on sgi is mpi
	   mpi_libs='-lmpi'

	   # set up sgi mpi defaults
	   if test -z "${MPI_LIB}" ; then
	       if test "${enable_32_bit}" = no ; then
		   MPI_LIB="${MPI_SGI}/usr/lib64"
	       else
		   MPI_LIB="${MPI_SGI}/usr/lib32"
	       fi
	   fi

       elif test "${with_mpi}" = mpich ; then

	   # no mpich support on SGI
	   AC_MSG_ERROR("We do not support mpich on the SGI yet!")

       fi

       # MPT (Message Passing Toolkit) for SGI vendor
       # implementation of MPI
       if test -z "${MPI_INC}" &&  test "${with_mpi}" = vendor ; then
	   MPI_INC="${MPT_SGI}/usr/include/"
       fi

       # add SGI MPT Specfic options
       if test "${with_mpi}" = vendor ; then
	   # define no C++ bindings
	   AC_DEFINE(MPI_NO_CPPBIND)
       fi

       #
       # end of communication package setup
       #

       #
       # setup lapack
       #

       if test "${with_lapack}" = vendor ; then
	   lapack_libs='-lcomplib.sgimath'
       fi

       #
       # end of lapack setup
       #

       #
       # gandolf, pcg and eospac requires -lfortran on the link line.
       #

       AC_MSG_CHECKING("libfortran requirements")
       if test -n "${vendor_gandolf}" || \
          test -n "${vendor_eospac}"  || \
          test -n "${vendor_pcg}" ; then
          LIBS="${LIBS} -lfortran"
          AC_MSG_RESULT("-lfortran added to LIBS")
       else
	   AC_MSG_RESULT("not needed")
       fi
       
       #
       # end of libfortran setup (gandolf, eospac, pcg)
       #

       #
       # pcg requires -lperfex on the link line.
       #

       AC_MSG_CHECKING("libperfex requirements")
       if test -n "${vendor_pcg}" ; then
          LIBS="${LIBS} -lperfex"
          AC_MSG_RESULT("-lperfex added to LIBS")
       else
	   AC_MSG_RESULT("not needed")
       fi
       
       #
       # end of libfortran setup (gandolf, eospac, pcg)
       #

       #
       # finalize vendors
       #
       AC_VENDOR_FINALIZE

       AC_DBS_SETUP_RPATH(-rpath, colon)

]) dnl irix

dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_DARWIN_ENVIRONMENT
dnl
dnl Configure draco build system Darwin-specific variables
dnl This function is called within AC_DBS_PLATFORM_ENVIRONMENT
dnl ***** NOT FULLY IMPLEMENTED *****
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_DARWIN_ENVIRONMENT], [dnl

       # dependency rules for IBM visual age compiler are complex
       if test "${with_cxx}" = ibm; then
	   DEPENDENCY_RULES='Makefile.dep.xlC.darwin'
       fi

       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       AC_DBS_SETUP_POSIX

       #   
       # LONG LONG on Linux
       #
       
       # always allow long long in strict ansi mode (if possible)
       
       if test -n "${STRICTFLAG}"; then

           case $CXX in

           # GNU g++
           g++) 
               AC_MSG_NOTICE([g++ -ansi option set to allow long long type!])
               STRICTFLAG="$STRICTFLAG -Wno-long-long"
               AC_MSG_NOTICE([g++ -ansi option set to allow long double type])
               STRICTFLAG="$STRICTFLAG -Wno-long-double"
           ;;
  	   ibm)	
	       AC_MSG_WARN("xlC set to allow long long")
	       STRICTFLAG="-qlanglvl=extended"
	       CFLAGS="${CFLAGS} -qlonglong"
	       CXXFLAGS="${CXXFLAGS} -qlonglong"
	   ;;
           # catchall
           *) 
               # do nothing
           ;;

           esac

       fi

       # 
       # end of LONG LONG setup
       #

       #
       # Setup communications packages
       #
       AC_DBS_SETUP_COMM(mpich)
	mpi_libs="-lmpich -lpmpich -lz"

       # 
       # setup lapack 
       #
       
       # we assume that the vendor option on linux is the install of
       # redhat rpms in /usr/lib; we don't worry about atlas because
       # that has already been defined

       if test "${with_lapack}" = vendor ; then
	   lapack_libs='-llapack -lblas'
       fi 

       # 
       # end of lapack setup
       # 

       # setup lf95 libs
       AC_DBS_LAHEY_ENVIRONMENT

       #
       # add libg2c to LIBS if lapack, gandolf, or pcg is used
       #
       AC_MSG_CHECKING("libg2c requirements")
       if test -n "${vendor_lapack}" || test -n "${vendor_pcg}" ||
	  test -n "${vendor_gandolf}"; then
	   
	   # Add g2c for various compilers
	   if test "${CXX}" = g++ ; then
	       LIBS="${LIBS} -lg2c"
	       AC_MSG_RESULT("-lg2c added to LIBS")
	   elif test "${CXX}" = icpc ; then
               AC_PATH_PROG(GCC_BIN, g++, null)
               GCC_BIN=`dirname ${GCC_BIN}`
               GCC_HOME=`dirname ${GCC_BIN}`
               GCC_LIB_DIR="${GCC_HOME}/lib"
	       LIBS="${LIBS} -L${GCC_LIB_DIR} -lg2c"
	       AC_MSG_RESULT("-lg2c added to LIBS")
           else
               AC_MSG_RESULT("not needed")
	   fi

       else
	   AC_MSG_RESULT("not needed")
       fi


       #
       # If dlopen is specified, 1) add libdl to LIBS; 
       # 2) add -fPIC to compile flags.
       #
       AC_MSG_CHECKING("libdl requirements")
       if test -n "${vendor_dlopen}" ; then
           if test "${enable_dlopen}" = yes ; then
               LIBS="${LIBS} -ldl"

               # if we are using g++ add fPIC
               if test "${CXX}" = g++ ; then
                   CFLAGS="${CFLAGS} -fPIC"
                   CXXFLAGS="${CXXFLAGS} -fPIC"
                   AC_MSG_RESULT("-ldl added to LIBS -fPIC added to compile flags")
               else
                   AC_MSG_RESULT("-ldl added to LIBS")
               fi

           else  
               AC_MSG_RESULT("not needed")
           fi

       else
           AC_MSG_RESULT("not needed")
       fi

       #
       # Set up fpe_trap for this platform.
       #
       AC_DEFINE(FPETRAP_DARWIN_PPC)

       #
       # finalize vendors
       #
       AC_VENDOR_FINALIZE

])

dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_SETUP_POSIX
dnl
dnl we do not do any posix source defines unless the user specifically
dnl requests them. 
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_SETUP_POSIX], [dnl

       if test "${with_posix:=no}" = yes ; then
	   with_posix='199309L'
       fi

       if test "${with_posix}" != no ; then
	   AC_DEFINE(_POSIX_SOURCE)
	   AC_DEFINE_UNQUOTED(_POSIX_C_SOURCE, $with_posix)
       fi

]) dnl setup_posix

dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_SETUP_COMM
dnl
dnl Setup communication packages
dnl
dnl default locations for mpi include/lib are:
dnl          /usr/local/mpich/include
dnl          /usr/local/mpich/lib
dnl to make life easy for CCS-2/4 users; needless to say,
dnl these can be overridden by --with-mpi-lib and --with-mpi-inc
dnl
dnl First argument is the default value for ${with_mpi} when this
dnl variable has the value 'vendor'.  
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_SETUP_COMM], [dnl

       # setup for mpi support, on linux vendor and mpich are one
       # and the same because there is no vendor for mpi on linux
        
       if test "${with_mpi}" = vendor ; then
	   with_mpi=$1
       fi

       # For CCS-2/4 users, we can also specify LAMPI in place of
       # mpich. 

       case ${with_mpi} in
       mpich)
	   # define mpi libs for mpich on linux
	   mpi_libs='-lmpich'
           ;;
       lampi | LAMPI | LA-MPI)
           with_mpi='LAMPI'
	   # define mpi libs for LAMPI on linux
	   mpi_libs='-lmpi'
           AC_MSG_CHECKING("mpirun -version")
           mpi_version=`mpirun -version`
           if test -n "`echo ${mpi_version} | grep -i LA-MPI`"; then
              AC_MSG_RESULT(${mpi_version})
           else
              AC_MSG_ERROR("Did not find LA-MPI version of mpirun.")
           fi
           ;;
       esac 
])


dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_SETUP_RPATH
dnl
dnl set rpath when building shared library executables
dnl
dnl We support two forms for RPATH support:
dnl 1) "-rpath dir1 -Xlinker -rpath dir2 ..."
dnl 2) "-rpath dir1:dir2:..."
dnl
dnl Some compilers/linkers use "R" instead of "rpath".  The option
dnl name is set from the 1st argument to this function.  The second
dnl argument specifies the list type as desribed above.
dnl
dnl $1 = rpath trigger.  One of "rpath" or "R"
dnl $2 = delimiter. One of "space", "nospace", or "colon"
dnl
dnl Some compilers require a pass-to-linker argument (ie. -Xlinker in
dnl g++).  These should be added to the rpath trigger argument, ie.
dnl
dnl AC_DBS_SETUP_RPATH('-Xlinker -rpath', space)
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_SETUP_RPATH], [dnl

       rptrigger=$1
       dilem=$2

       if test "${enable_shared}" = yes ; then

	   # turn off ranlib
	   RANLIB=':'

           if test "${dilem}" = "space"; then
	       RPATHA="${rptrigger} \${curdir}"
	       RPATHB="${rptrigger} \${curdir}/.."
	       RPATHC="${rptrigger} \${libdir}"
	       RPATH="${RPATHA} ${RPATHB} ${RPATHC} ${RPATH}"
           elif test "${dilem}" = "nospace"; then
	       RPATHA="${rptrigger}\${curdir}"
	       RPATHB="${rptrigger}\${curdir}/.."
	       RPATHC="${rptrigger}\${libdir}"
	       RPATH="${RPATHA} ${RPATHB} ${RPATHC} ${RPATH}"
           elif test "${dilem}" = "colon"; then
	       RPATH="${rptrigger} \${curdir}:\${curdir}/..:\${libdir} ${RPATH}"
           else
               AC_MSG_ERROR("Cannot determine what rpath format to use!")
	   fi
       fi

       # add vendors to rpath
       for vendor_dir in ${VENDOR_LIB_DIRS}; 
       do
           if test "${dilem}" = "space"; then
	       RPATH="${rptrigger} ${vendor_dir} ${RPATH}"
           elif test "${dilem}" = "nospace"; then
	       RPATH="${rptrigger}${vendor_dir} ${RPATH}"
           elif test "${dilem}" = "colon"; then
	       RPATH="${rptrigger} ${vendor_dir} ${RPATH}"
           else
               AC_MSG_ERROR("Cannot determine what rpath format to use!")
	   fi
       done

]) dnl setup_rpath

dnl-------------------------------------------------------------------------dnl
dnl end of ac_platforms.m4
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl ac_doxygen.m4
dnl
dnl Macros to help setup doxygen autodoc directories.
dnl
dnl Kelly Thompson
dnl 2004/03/30 16:41:22
dnl 1999/02/04 01:56:19
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl AC_SET_DEFAULT_OUTPUT
dnl-------------------------------------------------------------------------dnl
#
# Set the default location for doxygen output
#
AC_DEFUN([AC_SET_DEFAULT_OUTPUT], [dnl
   if test ${doxygen_output_top} = DEFAULT; then
       AC_SUBST(doxygen_output_top, "${prefix}/documentation")
   fi
])

dnl-------------------------------------------------------------------------dnl
dnl AC_AUTODOC_PACKAGE_TAGS
dnl
dnl  Collect tagfiles for pacakge-to-component dependencies
dnl-------------------------------------------------------------------------dnl
AC_DEFUN([AC_AUTODOC_PACKAGE_TAGS], [dnl

   # XXX Need to change COMPLINKS to generic doxygen list instead of
   # HTML for Latex compatability. Let doxygen insert the links
   AC_MSG_CHECKING([for documented sub-components of this package])
   COMP_LINKS=''
   TAGFILES=''
   DOXYGEN_TAGFILES=''
   components=''
   for item in `ls -1 ${package_top_srcdir}/src`; do
      if test -d ${package_top_srcdir}/src/${item}/autodoc; then
         dirname=`basename ${item}`
         components="${components} ${dirname}"
         COMP_LINKS="${COMP_LINKS} <li><a href=\"${dirname}/index.html\">${dirname}</a></li>"
         tagfile=${doxygen_output_top}/${dirname}.tag
         TAGFILES="${TAGFILES} ${tagfile}"
         DOXYGEN_TAGFILES="${DOXYGEN_TAGFILES} \"${tagfile} = ${dirname}\""
      fi
   done
   AC_MSG_RESULT(${components:-none})
   COMP_LINKS="<ul> $COMP_LINKS </ul>"

   # XXX TO DO: Add links to dependent packages on this page.
   PACKAGE_LINKS="<ul> </ul>"

   # Unique to package-level
   AC_SUBST(PACKAGE_LINKS)
   AC_SUBST(COMP_LINKS)

])


dnl-------------------------------------------------------------------------dnl
dnl AC_AUTODOC_COMPONENT_TAGS
dnl
dnl   Collect tagfiles for within-package component dependencies
dnl-------------------------------------------------------------------------dnl
#
# Build a list of tagfiles for other components of the same package
# and the _relative_ locations of the autodoc directories that they
# refer to.
#
# The relative path between component documentation in the same
# package is "../component" 
#
# These components are specified in AC_NEEDS_LIBS, and are stored
# in variable DEPENDENT_COMPONENTS. 
#
AC_DEFUN([AC_AUTODOC_COMPONENT_TAGS], [dnl

   components=''
   TAGFILES=''
   DOXYGEN_TAGFILES=''
   AC_MSG_CHECKING([for Doxygen component dependencies])
   for comp in ${DEPENDENT_COMPONENTS}; do
       components="${components} ${comp}"
       tagfile=${doxygen_output_top}/${comp}.tag
       DOXYGEN_TAGFILES="${DOXYGEN_TAGFILES} \"${tagfile} = ../${comp}\""
   done
   AC_MSG_RESULT([${components}])

])

dnl-------------------------------------------------------------------------dnl
dnl AC_AUTODOC_SUBST
dnl 
dnl   Do subsistutions on common AUTODOC variables
dnl-------------------------------------------------------------------------dnl
AC_DEFUN([AC_AUTODOC_SUBST], [dnl

   # Doxygen Input
   AC_SUBST(doxygen_input)
   AC_SUBST(doxygen_examples)

   # Doxygen Output
   AC_SUBST(doxygen_output_top)
   AC_SUBST(doxygen_html_output)
   AC_SUBST(doxygen_latex_output)

   # Other doxygen configuration
   AC_SUBST(DOXYGEN_TAGFILES)

   # For inclusion in header files and other html
   AC_SUBST(rel_package_html)

   # For makefiles for configuration:
   AC_SUBST(header_dir)
   AC_SUBST(autodoc_dir)

])

dnl-------------------------------------------------------------------------dnl
dnl AC_DRACO_AUTODOC
dnl
dnl  setup doxygen autodoc directories for COMPONENTS within a package
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DRACO_AUTODOC], [dnl

   # Get the default output location
   AC_SET_DEFAULT_OUTPUT

   # Define some package-level directories
   header_dir=${package_top_srcdir}/autodoc/html
   config_dir=${package_top_srcdir}/config

   abs_srcdir=`cd ${srcdir}; pwd`
   autodoc_dir=${abs_srcdir}/autodoc

   # For a component, the doxygen input is the srcdir and the examples
   # are in the tests
   AC_MSG_CHECKING([doxygen input directories])
   if test -d ${abs_srcdir}; then
      doxygen_input="${doxygen_input} ${abs_srcdir}"
   fi
   if test -d ${autodoc_dir}; then
      doxygen_input="${doxygen_input} ${autodoc_dir}"
   fi
   AC_MSG_RESULT(${doxygen_input})
   if test -d ${abs_srcdir}/test; then
      doxygen_examples=${abs_srcdir}/test
   fi

   # Set the package-level html output location
   package_html=${doxygen_output_top}/html

   # The local dir is different from the current dir.
   # localdir=`pwd`/autodoc

   # Set the component output locations.
   doxygen_html_output="${doxygen_output_top}/html/${package}"
   doxygen_latex_output="${doxygen_output_top}/latex/${package}"

   # Relative location of the package-level html output.
   adl_COMPUTE_RELATIVE_PATHS([doxygen_html_output:package_html:rel_package_html])

   # Get tags for other components in this package which this
   # component depends on
   AC_AUTODOC_COMPONENT_TAGS

   # find the release number
   number=$1
   AC_MSG_CHECKING("component release number")
   AC_MSG_RESULT($number)
   AC_SUBST(number)

   AC_AUTODOC_SUBST

   AC_CONFIG_FILES([autodoc/Makefile:${config_dir}/Makefile.autodoc.in \
                    autodoc/doxygen_config:${config_dir}/doxygen_config.in \
                    autodoc/header.html:${header_dir}/header.html.in \
                    autodoc/footer.html:${header_dir}/footer.html.in ])

])

dnl-------------------------------------------------------------------------dnl
dnl AC_PACKAGE_AUTODOC
dnl
dnl  setup doxygen autodoc directories for a PACKAGE
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_PACKAGE_AUTODOC], [dnl

   # Get the default output location
   AC_SET_DEFAULT_OUTPUT

   # Package-level directories
   header_dir=${srcdir}/html
   config_dir=${package_top_srcdir}/config

   abs_srcdir=`cd ${srcdir}; pwd`
   autodoc_dir=${abs_srcdir}

   # For the package, the input is the current directory, plus
   # configure/doc. There are no examples
   AC_MSG_CHECKING([for Doxygen input directories])
   doxygen_input="`pwd`"
   if test -d ${config_dir}/doc; then
      doxygen_input="${doxygen_input} ${config_dir}/doc"
   fi
   if test -d ${autodoc_dir}; then
      doxygen_input="${doxygen_input} ${autodoc_dir}"
   fi
   AC_MSG_RESULT(${doxygen_input})
   doxygen_examples=''

   # Component output locations
   doxygen_html_output="${doxygen_output_top}/html/"
   doxygen_latex_output="${doxygen_output_top}/latex/"

   # Relative location of the package-level html output.
   rel_package_html='.'

   AC_AUTODOC_PACKAGE_TAGS

   AC_AUTODOC_SUBST

   AC_CONFIG_FILES([doxygen_config:${config_dir}/doxygen_config.in])
   AC_CONFIG_FILES([Makefile:${config_dir}/Makefile.autodoc.in])
   AC_CONFIG_FILES([header.html:html/header.html.in])
   AC_CONFIG_FILES([footer.html:html/footer.html.in])
   AC_CONFIG_FILES([mainpage.dcc])

])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_doxygen.m4
dnl-------------------------------------------------------------------------dnl


dnl-------------------------------------------------------------------------dnl
dnl ac_utils.m4
dnl
dnl Macros to perform useful functions
dnl
dnl Mike Buksas
dnl-------------------------------------------------------------------------dnl

dnl Functions taken from:
dnl http://www.gnu.org/software/ac-archive/htmldoc/relpaths.html
dnl

AC_DEFUN([adl_COMPUTE_RELATIVE_PATHS],
[for _lcl_i in $1; do
  _lcl_from=\[$]`echo "[$]_lcl_i" | sed 's,:.*$,,'`
  _lcl_to=\[$]`echo "[$]_lcl_i" | sed 's,^[[^:]]*:,,' | sed 's,:[[^:]]*$,,'`
  _lcl_result_var=`echo "[$]_lcl_i" | sed 's,^.*:,,'`
  adl_RECURSIVE_EVAL([[$]_lcl_from], [_lcl_from])
  adl_RECURSIVE_EVAL([[$]_lcl_to], [_lcl_to])
  _lcl_notation="$_lcl_from$_lcl_to"
  adl_NORMALIZE_PATH([_lcl_from],['/'])
  adl_NORMALIZE_PATH([_lcl_to],['/'])
  adl_COMPUTE_RELATIVE_PATH([_lcl_from], [_lcl_to], [_lcl_result_tmp])
  adl_NORMALIZE_PATH([_lcl_result_tmp],["[$]_lcl_notation"])
  eval $_lcl_result_var='[$]_lcl_result_tmp'
done])


dnl adl_COMPUTE_RELATIVE_PATH(FROM, TO, RESULT)
dnl ===========================================
dnl Compute the relative path to go from $FROM to $TO and set the value
dnl of $RESULT to that value.  This function work on raw filenames
dnl (for instead it will considerate /usr//local and /usr/local as
dnl two distinct paths), you should really use adl_COMPUTE_REALTIVE_PATHS
dnl instead to have the paths sanitized automatically.
dnl
dnl For instance:
dnl    first_dir=/somewhere/on/my/disk/bin
dnl    second_dir=/somewhere/on/another/disk/share
dnl    adl_COMPUTE_RELATIVE_PATH(first_dir, second_dir, first_to_second)
dnl will set $first_to_second to '../../../another/disk/share'.
AC_DEFUN([adl_COMPUTE_RELATIVE_PATH],
[adl_COMPUTE_COMMON_PATH([$1], [$2], [_lcl_common_prefix])
adl_COMPUTE_BACK_PATH([$1], [_lcl_common_prefix], [_lcl_first_rel])
adl_COMPUTE_SUFFIX_PATH([$2], [_lcl_common_prefix], [_lcl_second_suffix])
$3="[$]_lcl_first_rel[$]_lcl_second_suffix"])

dnl adl_COMPUTE_COMMON_PATH(LEFT, RIGHT, RESULT)
dnl ============================================
dnl Compute the common path to $LEFT and $RIGHT and set the result to $RESULT.
dnl
dnl For instance:
dnl    first_path=/somewhere/on/my/disk/bin
dnl    second_path=/somewhere/on/another/disk/share
dnl    adl_COMPUTE_COMMON_PATH(first_path, second_path, common_path)
dnl will set $common_path to '/somewhere/on'.
AC_DEFUN([adl_COMPUTE_COMMON_PATH],
[$3=''
_lcl_second_prefix_match=''
while test "[$]_lcl_second_prefix_match" != 0; do
  _lcl_first_prefix=`expr "x[$]$1" : "x\([$]$3/*[[^/]]*\)"`
  _lcl_second_prefix_match=`expr "x[$]$2" : "x[$]_lcl_first_prefix"`
  if test "[$]_lcl_second_prefix_match" != 0; then
    if test "[$]_lcl_first_prefix" != "[$]$3"; then
      $3="[$]_lcl_first_prefix"
    else
      _lcl_second_prefix_match=0
    fi
  fi
done])

dnl adl_COMPUTE_SUFFIX_PATH(PATH, SUBPATH, RESULT)
dnl ==============================================
dnl Substrack $SUBPATH from $PATH, and set the resulting suffix
dnl (or the empty string if $SUBPATH is not a subpath of $PATH)
dnl to $RESULT.
dnl
dnl For instace:
dnl    first_path=/somewhere/on/my/disk/bin
dnl    second_path=/somewhere/on
dnl    adl_COMPUTE_SUFFIX_PATH(first_path, second_path, common_path)
dnl will set $common_path to '/my/disk/bin'.
AC_DEFUN([adl_COMPUTE_SUFFIX_PATH],
[$3=`expr "x[$]$1" : "x[$]$2/*\(.*\)"`])

dnl adl_COMPUTE_BACK_PATH(PATH, SUBPATH, RESULT)
dnl ============================================
dnl Compute the relative path to go from $PATH to $SUBPATH, knowing that
dnl $SUBPATH is a subpath of $PATH (any other words, only repeated '../'
dnl should be needed to move from $PATH to $SUBPATH) and set the value
dnl of $RESULT to that value.  If $SUBPATH is not a subpath of PATH,
dnl set $RESULT to the empty string.
dnl
dnl For instance:
dnl    first_path=/somewhere/on/my/disk/bin
dnl    second_path=/somewhere/on
dnl    adl_COMPUTE_BACK_PATH(first_path, second_path, back_path)
dnl will set $back_path to '../../../'.
AC_DEFUN([adl_COMPUTE_BACK_PATH],
[adl_COMPUTE_SUFFIX_PATH([$1], [$2], [_lcl_first_suffix])
$3=''
_lcl_tmp='xxx'
while test "[$]_lcl_tmp" != ''; do
  _lcl_tmp=`expr "x[$]_lcl_first_suffix" : "x[[^/]]*/*\(.*\)"`
  if test "[$]_lcl_first_suffix" != ''; then
     _lcl_first_suffix="[$]_lcl_tmp"
     $3="../[$]$3"
  fi
done])

dnl adl_RECURSIVE_EVAL(VALUE, RESULT)
dnl =================================
dnl Interpolate the VALUE in loop until it doesn't change,
dnl and set the result to $RESULT.
dnl WARNING: It's easy to get an infinite loop with some unsane input.
AC_DEFUN([adl_RECURSIVE_EVAL],
[_lcl_receval="$1"
$2=`(test "x$prefix" = xNONE && prefix="$ac_default_prefix"
     test "x$exec_prefix" = xNONE && exec_prefix="${prefix}"
     _lcl_receval_old=''
     while test "[$]_lcl_receval_old" != "[$]_lcl_receval"; do
       _lcl_receval_old="[$]_lcl_receval"
       eval _lcl_receval="\"[$]_lcl_receval\""
     done
     echo "[$]_lcl_receval")`])



dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/normpath.html
dnl
AC_DEFUN([adl_NORMALIZE_PATH],
[case ":[$]$1:" in
# change empty paths to '.'
  ::) $1='.' ;;
# strip trailing slashes
  :*[[\\/]]:) $1=`echo "[$]$1" | sed 's,[[\\/]]*[$],,'` ;;
  :*:) ;;
esac
# squeze repeated slashes
case ifelse($2,,"[$]$1",$2) in
# if the path contains any backslashes, turn slashes into backslashes
 *\\*) $1=`echo "[$]$1" | sed 's,\(.\)[[\\/]][[\\/]]*,\1\\\\,g'` ;;
# if the path contains slashes, also turn backslashes into slashes
 *) $1=`echo "[$]$1" | sed 's,\(.\)[[\\/]][[\\/]]*,\1/,g'` ;;
esac])




