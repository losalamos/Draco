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
dnl Service macros used in configure.in's throughout Draco.
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
dnl usage: configure.in
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
])

dnl-------------------------------------------------------------------------dnl
dnl AC_NEEDS_LIBS_TEST
dnl
dnl add DRACO-dependent libraries necessary for a package test
dnl usage: configure.in
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
dnl usage: in configure.in:
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
dnl usage: configure.in
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
dnl usage: configure.in
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
dnl usage: configure.in
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
dnl end of ac_conf.m4
dnl-------------------------------------------------------------------------dnl


dnl-------------------------------------------------------------------------dnl
dnl ac_dracoenv.m4
dnl
dnl Defines the Draco build system environment.  This is the main
dnl configure function.
dnl
dnl Thomas M. Evans
dnl 1999/02/04 01:56:21
dnl-------------------------------------------------------------------------dnl

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
   INSTALL_DATA='${INSTALL} -m 444'

   dnl
   dnl C4 OPERATIONS
   dnl

   # do the correct #defines
   if test "$with_c4" = scalar ; then
       AC_DEFINE(C4_SCALAR)
   elif test "$with_c4" = mpi ; then
       AC_DEFINE(C4_MPI)
   fi

   # if c4=mpi and with-mpi=no explicitly then 
   # define them (mpi gets set to vendor by default)
   if test "$with_c4" = mpi ; then
       if test "$with_mpi" = no ; then
	   with_mpi='vendor'
       fi
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

   # do draco standard headers

   AC_MSG_CHECKING("for draco standard headers")
   if test "${enable_draco_stdhdrs:=no}" != no ; then

       # if draco vendor is defined then use that include path
       if test -n "${DRACO_INC}" ; then
	   CPPFLAGS="${CPPFLAGS} -I${DRACO_INC}/stdheaders"

       # otherwise use the standard install location (includedir)
       else
	   CPPFLAGS="${CPPFLAGS} "'-I${includedir}/stdheaders'

       fi
       AC_MSG_RESULT("CPPFLAGS modified")

   else

       # we don't need draco stdheaders
       AC_MSG_RESULT("no") 

   fi

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
      [  --with-cxx[=gcc,sgi,kcc,compaq,guide]                                    
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
   dnl DRACO STANDARD HEADERS
   dnl

   dnl defines --enable-draco-stdhdrs
   AC_ARG_ENABLE(draco-stdhdrs,
      [  --enable-draco-stdhdrs  use draco standard headers (off by default)])

   dnl Doxygen options

   AC_ARG_ENABLE(latex-doc,
      [  --enable-latex-doc      built latex docs with doxygen (off by default)])

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

   # Check for working Lahey F90 compiler

   AC_CHECK_PROG(F90, lf95, lf95, none)
   if test "${F90}" = lf95 && ${F90} --version 2>&1 | grep "Lahey"
   then
       :
   else
       AC_MSG_ERROR([not found])
   fi
  
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
       F90FLAGS="--f95 --in --info --swm 2004,2006,2008,8202,8203,8204,8205,8206,8209,8220 ${F90FREE}"

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
	    F90FLAGS="-O5 -arch host -assume noaccuracy_sensitive -math_library fast -tune host ${F90FLAGS}"
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
   
   if test "${with_cxx}" = kcc ; then
       AC_CHECK_PROG(CXX, KCC, KCC)

       if test "${CXX}" = KCC ; then
	   CC='KCC --c'
	   AC_DRACO_KCC
       else
	   AC_MSG_ERROR("Did not find KCC compiler!")
       fi

   elif test "${with_cxx}" = sgi ; then
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

   elif test "${with_cxx}" = icc ; then 
       AC_CHECK_PROG(CXX, icc, icc)

       if test "${CXX}" = icc ; then
	   CC='icc'
	   AC_DRACO_INTEL_ICC
       else
	   AC_MSG_ERROR("Did not find Intel icc compiler!")
       fi

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

   # set the language to CPLUSPLUS
   AC_LANG_CPLUSPLUS

])

dnl-------------------------------------------------------------------------dnl
dnl KCC COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_KCC, [dnl

   AC_MSG_CHECKING("configuration of ${CXX}/${CC} compilers")

   # KCC SPECIFIC FLAGS
   dirstoclean='ti_files'

   # LINKER AND LIBRARY (AR)
   LD='${CXX}'
   AR='${CXX}'
   ARFLAGS='-o'
   ARLIBS='${DRACO_LIBS}'
   ARTESTLIBS='${PKG_LIBS} ${DRACO_TEST_LIBS} ${DRACO_LIBS}'

   # COMPILATION FLAGS

   # strict asci compliance
   if test "${enable_strict_ansi:=yes}" = yes ; then
       STRICTFLAG="--strict -D__KAI_STRICT"
   fi

   # --one_per flag
   if test "${enable_one_per:=yes}" = yes ; then
       ONEPERFLAG="--one_per"
   fi

   # optimization level
   if test "${enable_debug:=no}" = yes && \
      test "${with_opt:=0}" != 0 ; then
      CXXFLAGS="${CXXFLAGS} -g"
      CFLAGS="${CFLAGS} -g"
   fi
   CXXFLAGS="${CXXFLAGS} +K${with_opt:=0}"
   CFLAGS="${CFLAGS} +K${with_opt:=0}"

   # static linking option
   if test "${enable_static_ld}" = yes ; then
       LDFLAGS="${LDFLAGS} --static_libKCC -Bstatic"
   fi

   # Parallel build flag
   
   PARALLEL_FLAG="--parallel_build \${nj}"

   # final compiler additions
   CXXFLAGS="${CXXFLAGS} ${ONEPERFLAG}"

   AC_MSG_RESULT("KCC compiler flags set")
   
   dnl end of AC_DRACO_KCC
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
       GCC_LIB_DIR="${GCC_HOME}/lib"
   fi
   AC_MSG_RESULT("${GCC_LIB_DIR}")

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
       STRICTFLAG="-ansi -Wnon-virtual-dtor -Wreturn-type"
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
   # automatically); this, unfortunately, may become host dependent
   # in the future
   RPATH="${RPATH} -Xlinker -rpath ${GCC_LIB_DIR}"

   # static linking option
   if test "${enable_static_ld}" = yes ; then
       LDFLAGS="${LDFLAGS} -Bstatic"
   fi

   AC_MSG_RESULT("GNU g++ compiler flags set")

   dnl end of AC_DRACO_GNU_GCC
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
       AR="${CXX}"
       ARFLAGS="-shared -nocxxstd -expect_unresolved '*3td*' "
       ARFLAGS="${ARFLAGS} -expect_unresolved '*8_RWrwstd*' -o"
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

   # use implicit local template instantiation; this is the "GNU" like
   # option that puts manually instantiated templates in the 
   # repository with external linkage and automatic templates in 
   # the object file with internal linkage
   CXXFLAGS="${CXXFLAGS} -timplicit_local"

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
dnl Intel icc COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_INTEL_ICC, [dnl

   AC_MSG_CHECKING("configuration of ${CXX}/${CC} compilers")

   # icc SPECIFIC FLAGS

   # LINKER AND LIBRARY
   LD='${CXX}'

   # if shared then ar is icc
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
	   icc_opt_flags="-g -O${with_opt} -Ob1 -ip"
       else
	   icc_opt_flags="-O${with_opt} -Ob1"
       fi

   #set up compiler when not optimized (turn off inlining with -Ob0)
   else

       # turn on debug by default
       if test "${enable_debug:=yes}" = yes ; then
	   icc_opt_flags="-g -O0 -Ob0"
       else
	   icc_opt_flags="-O0 -Ob0"
       fi

   fi
   
   # set the cxx and c flags
   CXXFLAGS="${CXXFLAGS} ${icc_opt_flags}"
   CFLAGS="${CFLAGS} ${icc_opt_flags}"

   # static linking option
   if test "${enable_static_ld}" = yes ; then
       LDFLAGS="${LDFLAGS} -static"
   fi

   AC_MSG_RESULT("icc compiler flags set")
   
   dnl end of AC_DRACO_INTEL_ICC
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

       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       # we do not do any posix source defines unless the user
       # specifically requests them
       if test "${with_posix:=no}" = yes ; then
	   with_posix='199309L'
       fi

       if test "${with_posix}" != no ; then
	   AC_DEFINE(_POSIX_SOURCE)
	   AC_DEFINE_UNQUOTED(_POSIX_C_SOURCE, $with_posix)
       fi

       #
       # setup linux strict if the compiler is KCC (also turn off the
       # warnings about long long being non-standard)
       #
       if test "${CXX}" = KCC && test -n "${STRICTFLAG}" ; then
	   AC_MSG_WARN("Linux KCC strict option set to allow long long type!")
	   STRICTFLAG="--linux_strict -D__KAI_STRICT --diag_suppress 450"
       fi

       #
       # add thread safety if we are using KCC on linux
       #
       if test "${CXX}" = KCC ; then
	   CFLAGS="--thread_safe ${CFLAGS}"
	   CXXFLAGS="--thread_safe ${CXXFLAGS}"
	   ARFLAGS="--thread_safe ${ARFLAGS}"
	   LDFLAGS="--thread_safe ${LDFLAGS}"
       fi

       #
       # setup communication packages
       #
       
       # the default locations for mpi include/lib are:
       #   /usr/local/mpich/include
       #   /usr/local/mpich/lib
       # to make life easy for CCS-2/4 users; needless to say,
       # these can be overridden by --with-mpi-lib and --with-mpi-inc

       # setup for mpi support, on linux vendor and mich are one
       # and the same because there is no vendor for mpi on linux
        
       if test "${with_mpi}" = vendor ; then
	   with_mpi='mpich'
       fi

       if test "${with_mpi}" = mpich ; then

	   # define mpi libs for mpich on linux
	   mpi_libs='-lmpich'

	   # if /usr/local/mpich/lib exists use it by default;
	   # this is set as the default for the CCS-2/4 network;
	   # it may not be appropriate on other LINUX networks;
	   # in those cases, override with --with-mpi-lib
	   if test -z "${MPI_LIB}" && test -d "/usr/local/mpich/lib"; then
	       MPI_LIB='/usr/local/mpich/lib'
	   fi

	   # set the default include location on LINUX to
	   # /usr/local/mpich/include; this is specific to the CCS-2/4
	   # LINUX network; to override on other systems use
	   # --with-mpi-inc on the configure line

	   # if MPI_INC is undefined then define it
	   if test -z "${MPI_INC}" && test -d "/usr/local/mpich/include"; then
	       MPI_INC='/usr/local/mpich/include'
	   fi

       fi

       #
       # end of communication package setup
       #  

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

       #
       # setup eospac
       #
       
       AC_MSG_CHECKING("for extra eospac library requirements.")
       if test -n "${vendor_eospac}"; then
	   extra_eospac_libs="-L/usr/local/lf9562/lib -lfj9i6 -lfj9e6 -lfj9f6 -lfst -lfccx86_6a"
           LIBS="${LIBS} ${extra_eospac_libs}"
           AC_MSG_RESULT("${extra_eospac_libs}")
       else
           AC_MSG_RESULT("none.")
       fi

       #
       # end of eospac
       #

       #
       # add libg2c to LIBS if lapack, gandolf, or pcg is used
       #
       AC_MSG_CHECKING("libg2c requirements")
       if test -n "${vendor_lapack}" || test -n "${vendor_pcg}" ||
	  test -n "${vendor_gandolf}"; then
	   
	   # Add g2c for various compilers
	   if test "${CXX}" = KCC ; then
	       LIBS="${LIBS} --backend -lg2c"
	       AC_MSG_RESULT("--backend -lg2c added to LIBS")
	   elif test "${CXX}" = g++ ; then
	       LIBS="${LIBS} -lg2c"
	       AC_MSG_RESULT("-lg2c added to LIBS")
	   elif test "${CXX}" = icc ; then
               AC_PATH_PROG(GCC_BIN, g++, null)
               GCC_BIN=`dirname ${GCC_BIN}`
               GCC_HOME=`dirname ${GCC_BIN}`
               GCC_LIB_DIR="${GCC_HOME}/lib"
	       LIBS="${LIBS} -L${GCC_LIB_DIR} -lg2c"
	       AC_MSG_RESULT("-lg2c added to LIBS")
	   fi

       else
	   AC_MSG_RESULT("not needed")
       fi

       #
       # finalize vendors
       #
       AC_VENDOR_FINALIZE

       # set rpath when building shared library executables
       if test "${enable_shared}" = yes ; then

	   # turn off ranlib
	   RANLIB=':'

	   # the g++/icc rpath needs Xlinker in front of it
	   if test "${CXX}" = g++ || test "${CXX}" = icc; then
	       RPATHA="-Xlinker -rpath \${curdir}"
	       RPATHB="-Xlinker -rpath \${curdir}/.."
	       RPATHC="-Xlinker -rpath \${libdir}"
	       RPATH="${RPATHA} ${RPATHB} ${RPATHC} ${RPATH}"
	   else
	       RPATH="-rpath \${curdir}:\${curdir}/..:\${libdir} ${RPATH}"
	   fi

       fi

       # add vendors to rpath
       for vendor_dir in ${VENDOR_LIB_DIRS}; 
       do
	   # if we are using gcc/icc then add xlinker
	   if test "${CXX}" = g++ || test "${CXX}" = icc; then
	       RPATH="-Xlinker -rpath ${vendor_dir} ${RPATH}"

	   # else we just add the rpath
	   else
	       RPATH="-rpath ${vendor_dir} ${RPATH}"
	   fi
       done

       # add the intel math library for better performance when
       # compiling with intel
       if test "${CXX}" = icc; then
	   LIBS="$LIBS -limf"
       fi
   ;;

   # ***********
   # CYGWIN SETUP
   # ***********
   i686-pc-cygwin)

       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       # we do not do any posix source defines unless the user
       # specifically requests them
       if test "${with_posix:=no}" = yes ; then
	   with_posix='199309L'
       fi

       if test "${with_posix}" != no ; then
	   AC_DEFINE(_POSIX_SOURCE)
	   AC_DEFINE_UNQUOTED(_POSIX_C_SOURCE, $with_posix)
       fi

       #
       # setup linux strict if the compiler is KCC (also turn off the
       # warnings about long long being non-standard)
       #
       if test "${CXX}" = KCC && test -n "${STRICTFLAG}" ; then
	   AC_MSG_WARN("Linux KCC strict option set to allow long long type!")
	   STRICTFLAG="--linux_strict -D__KAI_STRICT --diag_suppress 450"
       fi

       #
       # add thread safety if we are using KCC on linux
       #
       if test "${CXX}" = KCC ; then
	   CFLAGS="--thread_safe ${CFLAGS}"
	   CXXFLAGS="--thread_safe ${CXXFLAGS}"
	   ARFLAGS="--thread_safe ${ARFLAGS}"
	   LDFLAGS="--thread_safe ${LDFLAGS}"
       fi

       #
       # setup communication packages
       #
       
       # the default locations for mpi include/lib are:
       #   /usr/local/mpich/include
       #   /us	   # define mpi libs for mpich on linux
	   mpi_libs='-lmpich'

	   dnl ifr/local/mpich/lib
       dnl to make life easy for CCS-2/4 users; needless to say,
       dnl these can be overridden by --with-mpi-lib and --with-mpi-inc

       dnl setup for mpi support, on linux vendor and mich are one
       dnl and the same because there is no vendor for mpi on linux
        
       if test "${with_mpi}" = vendor ; then
	   with_mpi='mpich'
       fi

       if test "${with_mpi}" = mpich ; then

	   dnl define mpi libs for mpich on linux
	   mpi_libs='-lmpich'

	   dnl if /usr/local/mpich/lib exists use it by default;
	   dnl this is set as the default for the CCS-2/4 network;
	   dnl it may not be appropriate on other LINUX networks;
	   dnl in those cases, override with --with-mpi-lib
	   if test -z "${MPI_LIB}" && test -d "/usr/local/mpich/lib"; then
	       MPI_LIB='/usr/local/mpich/lib'
	   fi

	   dnl set the default include location on LINUX to
	   dnl /usr/local/mpich/include; this is specific to the CCS-2/4
	   dnl LINUX network; to override on other systems use
	   dnl --with-mpi-inc on the configure line

	   dnl if MPI_INC is undefined then define it
	   if test -z "${MPI_INC}" && test -d "/usr/local/mpich/include"; then
	       MPI_INC='/usr/local/mpich/include'
	   fi

       fi

       dnl
       dnl end of communication package setup
       dnl  

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

       dnl
       dnl setup eospac
       dnl
       
       AC_MSG_CHECKING("for extra eospac library requirements.")
       if test -n "${vendor_eospac}"; then
	   extra_eospac_libs="-L/usr/local/lf9562/lib -lfj9i6 -lfj9e6 -lfj9f6 -lfst -lfccx86_6a"
           LIBS="${LIBS} ${extra_eospac_libs}"
           AC_MSG_RESULT("${extra_eospac_libs}")
       else
           AC_MSG_RESULT("none.")
       fi

       dnl
       dnl end of eospac
       dnl

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

       # set rpath when building shared library executables
       if test "${enable_shared}" = yes ; then

	   # turn off ranlib
	   RANLIB=':'

	   # the g++/icc rpath needs Xlinker in front of it
	   if test "${CXX}" = g++ || test "${CXX}" = icc; then
	       RPATHA="-Xlinker -rpath \${curdir}"
	       RPATHB="-Xlinker -rpath \${curdir}/.."
	       RPATHC="-Xlinker -rpath \${libdir}"
	       RPATH="${RPATHA} ${RPATHB} ${RPATHC} ${RPATH}"
	   else
	       RPATH="-rpath \${curdir}:\${curdir}/..:\${libdir} ${RPATH}"
	   fi

       fi

       # add vendors to rpath
       for vendor_dir in ${VENDOR_LIB_DIRS}; 
       do
	   # if we are using gcc/icc then add xlinker
	   if test "${CXX}" = g++ || test "${CXX}" = icc; then
	       RPATH="-Xlinker -rpath ${vendor_dir} ${RPATH}"

	   # else we just add the rpath
	   else
	       RPATH="-rpath ${vendor_dir} ${RPATH}"
	   fi
       done

       # add the intel math library for better performance when
       # compiling with intel
       if test "${CXX}" = icc; then
	   LIBS="$LIBS -limf"
       fi
   ;;

   # *********
   # SGI SETUP
   # *********
   mips-sgi-irix6.*)
   
       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       # posix source defines, by default we set posix on 
       if test "${with_posix:=yes}" = yes ; then
	   with_posix='199309L'
       fi

       if test "${with_posix}" != no ; then
	   AC_DEFINE_UNQUOTED(_POSIX_C_SOURCE, $with_posix)
	   AC_DEFINE(_POSIX_SOURCE)
       fi

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
       # adjust strict flag for KCC
       #

       if test "${CXX}" = KCC && test "${enable_strict_ansi}" = yes ; then
	   
	   # if integer type is long long or vendor mpi is on then we
	   # need to allow long long type
	   if test "${with_mpi}" = vendor; then 
		   AC_MSG_WARN("KCC strict option set to allow long long type")
		   STRICTFLAG="${STRICTFLAG} --diag_suppress 450"
	   fi

       fi

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
       # gandolf and eospac requires -lfortran on the link line.
       #

       AC_MSG_CHECKING("libfortran requirements")
       if test -n "${vendor_gandolf}" || test -n "${vendor_eospac}" ; then
          LIBS="${LIBS} -lfortran"
          AC_MSG_RESULT("-lfortran added to LIBS")
       else
	   AC_MSG_RESULT("not needed")
       fi
       
       #
       # end of gandolf/libfortran setup
       #

       #
       # finalize vendors
       #
       AC_VENDOR_FINALIZE

       # set rpath when building shared library executables
       if test "${enable_shared}" = yes; then

	   # the g++ rpath needs Xlinker in front of it
	   if test "${CXX}" = g++; then
	       RPATHA="-Xlinker -rpath \${curdir}"
	       RPATHB="-Xlinker -rpath \${curdir}/.."
	       RPATHC="-Xlinker -rpath \${libdir}"
	       RPATH="${RPATHA} ${RPATHB} ${RPATHC} ${RPATH}"
	   else
	       RPATH="-rpath \${curdir}:\${curdir}/..:\${libdir} ${RPATH}"
	   fi

       fi

       # add vendors to rpath
       for vendor_dir in ${VENDOR_LIB_DIRS}; 
       do
	   # if we are using gcc then add xlinker
	   if test "${CXX}" = g++; then
	       RPATH="-Xlinker -rpath ${vendor_dir} ${RPATH}"

	   # else we just add the rpath
	   else
	       RPATH="-rpath ${vendor_dir} ${RPATH}"
	   fi
       done
   ;;

   # ******************
   # TRU64 COMPAQ SETUP
   # ******************
   alpha*-dec-osf*)
   
       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       # posix source defines, by default we set posix off
       if test "${with_posix:=no}" = yes ; then
	   with_posix='199309L'
       fi

       if test "${with_posix}" != no ; then
	   AC_DEFINE_UNQUOTED(_POSIX_C_SOURCE, $with_posix)
	   AC_DEFINE(_POSIX_SOURCE)
       fi

       #
       # setup communication packages
       #

       # setup vendor mpi
       if test "${with_mpi}" = vendor ; then

	   # define mpi libraries
	   mpi_libs='-lmpi'
       
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

       if test "${with_lapack}" = vendor ; then
	   lapack_libs='-lcxmlp -lcxml'
       fi

       #
       # end of lapack setup
       #

       #
       # gandolf, eospac, pcg require -lfor on the link line.
       #

       AC_MSG_CHECKING("libfortran requirements")
       if test -n "${vendor_gandolf}" || test -n "${vendor_eospac}" ||
          test -n "${vendor_pcg}"; then
          LIBS="${LIBS} -lfor"
          AC_MSG_RESULT("-lfor added to LIBS")
       else
	   AC_MSG_RESULT("not needed")
       fi

       #
       # end of gandolf/libfortran setup
       #

       #
       # finalize vendors
       #
       AC_VENDOR_FINALIZE

       # set rpath when building shared library executables
       if test "${enable_shared}" = yes; then

	   # turn off ranlib
	   RANLIB=':'

	   # the g++ rpath needs Xlinker in front of it
	   if test "${CXX}" = g++; then
	       RPATHA="-Xlinker -rpath \${curdir}"
	       RPATHB="-Xlinker -rpath \${curdir}/.."
	       RPATHC="-Xlinker -rpath \${libdir}"
	       RPATH="${RPATHA} ${RPATHB} ${RPATHC} ${RPATH}"
	   else
	       RPATH="-rpath \${curdir}:\${curdir}/..:\${libdir} ${RPATH}"
	   fi

       fi

       # add vendors to rpath
       for vendor_dir in ${VENDOR_LIB_DIRS}; 
       do
	   # if we are using gcc then add xlinker
	   if test "${CXX}" = g++; then
	       RPATH="-Xlinker -rpath ${vendor_dir} ${RPATH}"

	   # else we just add the rpath
	   else
	       RPATH="-rpath ${vendor_dir} ${RPATH}"
	   fi
       done
   ;;

   # *************
   # IBM AIX SETUP
   # *************
   *ibm-aix*)

       # dependency rules for IBM visual age compiler are complex
       if test "${with_cxx}" = asciwhite || test "${with_cxx}" = ibm; then
	   DEPENDENCY_RULES='Makefile.dep.xlC'
       fi
   
       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       # posix source defines, by default we set posix off
       if test "${with_posix:=no}" = yes ; then
	   with_posix='199309L'
       fi

       if test "${with_posix}" != no ; then
	   AC_DEFINE_UNQUOTED(_POSIX_C_SOURCE, $with_posix)
	   AC_DEFINE(_POSIX_SOURCE)
       fi

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

       # RPATH is derived from -L, don't need explicit setup

       # do shared specific stuff
       if test "${enable_shared}" = yes ; then
	   # turn off ranlib
	   RANLIB=':'
       fi
   ;;

   # *****************
   # SUN/SOLARIS SETUP
   # *****************
   sparc-sun-solaris2.*)
   
       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       # posix source defines, by default we set poaix on 
       if test "${with_posix:=yes}" = yes ; then
	   with_posix='199309L'
       fi

       if test "${with_posix}" != no ; then
	   AC_DEFINE_UNQUOTED(_POSIX_C_SOURCE, $with_posix)
	   AC_DEFINE(_POSIX_SOURCE)
       fi

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

       # set -R when building shared library executables
       if test "${enable_shared}" = yes; then
	   RPATH="-R \${curdir}:\${curdir}/..:\${libdir} ${RPATH}"
	   RANLIB=':'
       fi

       # add vendors to rpath
       for vendor_dir in ${VENDOR_LIB_DIRS}; 
       do
	   # if we are using gcc/icc then add xlinker
	   if test "${CXX}" = g++ || test "${CXX}" = icc; then
	       RPATH="-Xlinker -R ${vendor_dir} ${RPATH}"

	   # else we just add the rpath
	   else
	       RPATH="-R ${vendor_dir} ${RPATH}"
	   fi
       done
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
dnl end of ac_platforms.m4
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
      [  --with-mpi=[vendor,mpich] 
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
      [  --with-aztec=[lib]      determine the aztec lib (aztec is the default])
 
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
      [  --with-gsl=[lib]      determine the gsl lib (gsl is the default])
 
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
	   AC_VENDORLIB_SETUP(vendor_gsl, -L${GSL_LIB} -l${with_gsl})
       elif test -z "${GSL_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gsl, -l${with_gsl})
       fi

       # add GSL directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${GSL_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${GSL_INC}"

   fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_GSLCBLAS_SETUP
dnl
dnl GSLCBLAS SETUP (on by default)
dnl GSLCBLAS is a required vendor
dnl
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_GSLCBLAS_SETUP], [dnl

   dnl define --with-gslcblas
   AC_ARG_WITH(gslcblas,
      [  --with-gslcblas=[lib]      determine the gslcblas lib (gslcblas is the default])
 
   dnl define --with-gslcblas-inc
   AC_WITH_DIR(gslcblas-inc, GSLCBLAS_INC, \${GSLCBLAS_INC_DIR},
	       [tell where GSLCBLAS includes are])

   dnl define --with-gslcblas-lib
   AC_WITH_DIR(gslcblas-lib, GSLCBLAS_LIB, \${GSLCBLAS_LIB_DIR},
	       [tell where GSLCBLAS libraries are])

   # set default value of gslcblas includes and libs
   if test "${with_gslcblas:=gslcblas}" = yes ; then
       with_gslcblas='gslcblas'
   fi

   # determine if this package is needed for testing or for the 
   # package
   vendor_gslcblas=$1

])


AC_DEFUN([AC_GSLCBLAS_FINALIZE], [dnl

   # set up the libraries and include path
   if test "${vendor_gslcblas}"; then

       # include path
       if test -n "${GSLCBLAS_INC}"; then 
	   # add to include path
	   VENDOR_INC="${VENDOR_INC} -I${GSLCBLAS_INC}"
       fi

       # library path
       if test -n "${GSLCBLAS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gslcblas, -L${GSLCBLAS_LIB} -l${with_gslcblas})
       elif test -z "${GSLCBLAS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_gslcblas, -l${with_gslcblas})
       fi

       # add GSLCBLAS directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${GSLCBLAS_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${GSLCBLAS_INC}"

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
      [  --with-trilinos=[lib]    determine the trilinos implementation (aztecoo is default])
 
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
	   VENDOR_INC="${VENDOR_INC} -I${TRILINOS_INC}amesos    -I${TRILINOS_INC}aztecoo  -I${TRILINOS_INC}epetra"
	   VENDOR_INC="${VENDOR_INC} -I${TRILINOS_INC}epetraext -I${TRILINOS_INC}ifpack   -I${TRILINOS_INC}komplex"
	   VENDOR_INC="${VENDOR_INC} -I${TRILINOS_INC}nox       -I${TRILINOS_INC}triutils -I${TRILINOS_INC}y12m"
       fi

       # library path
       if test -n "${TRILINOS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_trilinos, -L${TRILINOS_LIB} -l${with_trilinos} -lepetra -ltriutils -ly12m)
       elif test -z "${TRILINOS_LIB}" ; then
	   AC_VENDORLIB_SETUP(vendor_trilinos, -l${with_trilinos} -lepetra -ltriutils -ly12m)
       fi

       # add TRILINOS directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${TRILINOS_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${TRILINOS_INC}"

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

   AC_AZTEC_FINALIZE
   AC_PCG_FINALIZE

   AC_LAPACK_FINALIZE
   AC_EOSPAC5_FINALIZE
   AC_GANDOLF_FINALIZE
   AC_SPRNG_FINALIZE
   AC_GRACE_FINALIZE
   AC_METIS_FINALIZE

   AC_GSL_FINALIZE
   AC_GSLCBLAS_FINALIZE

   AC_MPI_FINALIZE

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
   AC_GSLCBLAS_SETUP(pkg)
   AC_TRILINOS_SETUP(pkg)
   AC_METIS_SETUP(pkg)
   AC_LAPACK_SETUP(pkg)
   AC_GANDOLF_SETUP(pkg)
   AC_EOSPAC5_SETUP(pkg)
   AC_GRACE_SETUP(pkg)
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_vendors.m4
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
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_local.m4
dnl-------------------------------------------------------------------------dnl

