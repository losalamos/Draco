dnl-------------------------------------------------------------------------dnl
dnl File  : draco/config ac_dracoenv.m4
dnl Author: Thomas M. Evans
dnl Date  : 1999/02/04 01:56:21
dnl
dnl Defines the Draco build system environment.  This is the main
dnl configure function.
dnl
dnl-------------------------------------------------------------------------dnl
dnl $Id$
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
      dnl We need to add the rpath to the LDFLAGS instead of RPATH
      dnl because configure fails AC_CHECK_SIZEOF if sltport is not
      dnl available.
      AC_MSG_CHECKING("for LDFLAGS mods for stlport")
      LDFLAGS="${stlport_xlinker} -rpath ${with_stlport}/lib ${LDFLAGS}"
      AC_MSG_RESULT("Added ${stlport_xlinker} -rpath ${with_stlport}/lib to LDFLAGS")

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

