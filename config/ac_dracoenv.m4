dnl-------------------------------------------------------------------------dnl
dnl ac_dracoenv.m4
dnl puts together the DRACO environments given the arguments from 
dnl vendors (ac_vendors.m4) and DRACO (ac_dracoarg.m4)
dnl
dnl $Id$
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl AC_DRACO_ENV
dnl
dnl usage: configure.in
dnl puts together the DRACO compile-time environments
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
   dnl RUN AC_PROG_INSTALL
   dnl

   AC_PROG_INSTALL

   dnl
   dnl NMLGEN PATH
   dnl

   NMLGEN=${libexecdir}/nmlgen

   dnl
   dnl C4 OPERATIONS
   dnl

   # do the correct #defines
   if test "$with_c4" = scalar ; then
       AC_DEFINE(C4_SCALAR)
   elif test "$with_c4" = mpi ; then
       AC_DEFINE(C4_MPI)
   elif test "$with_c4" = shmem ; then
       AC_DEFINE(C4_SHMEM)
   fi

   # if c4=mpi or shmem and with-mpi and enable-shmem are
   # set to no explicitly then define them (mpi gets set to   
   # vendor by default)
   if test "$with_c4" = mpi ; then
       if test "$with_mpi" = no ; then
	   with_mpi='vendor'
       fi
   elif test "$with_c4" = shmem ; then
       if test "$enable_shmem" = no ; then
	   enable_shmem='yes'
       fi
   fi

   # now set up the platform-independent comm directories
   AC_COMM_SET
   
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

   # set the POSIX source level
   if test "${with_posix}" = yes ; then
       AC_DEFINE(_POSIX_C_SOURCE, "199309L")
       AC_DEFINE(_POSIX_SOURCE)
   elif test "${with_posix:=199309L}" != no ; then
       AC_DEFINE_UNQUOTED(_POSIX_C_SOURCE, $with_posix)
       AC_DEFINE(_POSIX_SOURCE)
   fi
   
   dnl
   dnl TOOL CHECKS
   dnl
   
   dnl check for and assign the path to python
   AC_PATH_PROG(PYTHON_PATH, python, null)
   if test "${PYTHON_PATH}" = null ; then
       AC_MSG_ERROR("No valid Python found!")
   fi

   dnl
   dnl COMPILER SETUPS
   dnl

   dnl first find the host
   AC_CANONICAL_HOST

   dnl determine which compiler we are using

   # do tests of --with-cxx, see if the compiler exists and then call
   # the proper setup function, we default to the GNU EGCS compiler 
   # that is defined by g++ (c++) and gcc under EGCS
   
   if test "${with_cxx}" = kcc ; then
       AC_CHECK_PROG(CXX, KCC, KCC)
       if test "${CXX}" = KCC ; then
	   CC='KCC --c'
	   AC_DRACO_KCC
       else
	   AC_PROG_CXX
	   AC_PROG_CC
	   if test "${CXX}" = CC && test "${CC}" = cc ; then
	       AC_DRACO_CC
	   else
	       AC_DRACO_EGCS
	   fi
       fi
   elif test "${with_cxx}" = cc ; then
       AC_CHECK_PROG(CXX, CC, CC)
       AC_CHECK_PROG(CC, cc, cc)
       if test "${CXX}" = CC && test "${CC}" = cc ; then
	   AC_DRACO_CC
       else 
	   AC_PROG_CXX
	   AC_PROG_CC
	   AC_DRACO_EGCS
       fi
   elif test "${with_cxx}" = egcs ; then
       AC_PROG_CXX
       AC_PROG_CC
       AC_DRACO_EGCS
   fi
   
   # check to see that we have a C++ compiler defined, throw an error
   # if not
   if test "${CXX}" = KCC || test "${CXX}" = CC || \
      test "${CXX}" = g++ || test "${CXX}" = c++ ; then
	   found_cxx='good'
   fi
   
   if test "${found_cxx}" != good ; then
       AC_MSG_ERROR("No valid C++ Compiler Found!")
   fi
	
   dnl check for ranlib
   AC_PROG_RANLIB

   dnl setup the system-specific stuff

   # systems setup
   case $host in
   Linux)
       LD='${CXX}'
       AR="ar"
       ARFLAGS="cr"
       ARLIBS=""
   ;;
   mips-sgi-irix6.*)
       # RANLIB TAG ON SGI
       RANLIB=':'

       # BIT COMPILER FLAGS ON SGI
       if test "${enable_32_bit:=no}" = yes ; then
	   CXXFLAGS="-n32 ${CXXFLAGS}"
	   CFLAGS="-n32 ${CFLAGS}"
	   LDFLAGS="-n32 ${LDFLAGS}"
       else 
	   CXXFLAGS="-64 ${CXXFLAGS}"
	   CFLAGS="-64 ${CFLAGS}"
	   LDFLAGS="-64 ${LDFLAGS}"
       fi

       # MIPS INSTRUCTIONS ON SGI
       # this is different depending upon the compiler
       if test "${with_cxx}" = kcc ; then
	   CXXFLAGS="-mips${with_mips:=4} --backend -r10000 ${CXXFLAGS}"
	   CFLAGS="-mips${with_mips:=4} -r10000 ${CFLAGS}"
	   LDFLAGS="-mips${with_mips:=4} ${LDFLAGS}"
       elif test "${with_cxx}" = cc ; then
	   CXXFLAGS="-mips${with_mips:=4} -r10000 ${CXXFLAGS}"
	   CFLAGS="-mips${with_mips:=4} -r10000 ${CFLAGS}"
	   LDFLAGS="-mips${with_mips:=4} ${LDFLAGS}"
       fi

       # MPT (Message Passing Toolkit) for SGI vendor
       # implementation of MPI and SHMEM
       if test -z "${MPI_INC}" &&  test "${with_mpi}" = vendor ; then
	   MPI_INC="${MPT_SGI}/usr/include/"
	   AC_DEFINE_UNQUOTED(MPI_INC, ${MPI_INC})	   
	   MPI_H="\"${MPI_INC}mpi.h\""
	   AC_DEFINE_UNQUOTED(MPI_H, ${MPI_H})
       fi

       if test -z "${SHMEM_INC}" && test "${enable_shmem}" = yes ; then
	   SHMEM_INC="${MPT_SGI}/usr/include/mpp/"
	   AC_DEFINE_UNQUOTED(SHMEM_INC, ${SHMEM_INC})
	   SHMEM_H="\"${SHMEM_INC}shmem.h\""
	   AC_DEFINE_UNQUOTED(SHMEM_H, ${SHMEM_H})
       fi

       #PCGLIB EXTRAS
       if test "${enable_pcglib}" = yes ; then
	   AC_VENDORLIB_SETUP(vendor_pcglib, -lcomplib.sgimath)
       fi

       # set rpath when building shared library executables
       if test "${enable_shared}" = yes; then
	   DRACO_LIBS="-rpath \${libdir} ${DRACO_LIBS}"
	   PKG_LIBS="-rpath .. ${PKG_LIBS}"
       fi
   ;;
   sparc-sun-solaris2.*)
       # MPICH LIBRARY EXTRAS
       if test "${with_mpi}" = mpich ; then
	   AC_VENDORLIB_SETUP(vendor_mpi, -lsocket -lnsl)
       fi

       # PCGLIB EXTRAS
       if test "${enable_pcglib}" = yes ; then
	   pcg_sun_libs='-llapack -lblas -lF77 -lM77 -lsunmath'
	   AC_VENDORLIB_SETUP(vendor_pcglib, ${pcg_sun_libs})
       fi

       # set -R when building shared library executables
       if test "${enable_shared}" = yes; then
	   DRACO_LIBS="-R \${libdir} ${DRACO_LIBS}"
	   PKG_LIBS="-R .. ${PKG_LIBS}"
	   RANLIB=':'
       fi
   ;;
   *)
   esac

   # add system specific libraries
   LIBS='-lm'

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

   # if this is a parallel build, mark it the tests scalar
   if test "${with_c4}" = scalar ; then
       test_scalar="scalar"
   fi

   # define the TESTFLAGS, for parallel runs the processor will be
   # added later in the Makefile

   if test "${test_scalar}" = scalar ; then
       test_flags="--${test_exe:=binary}"
   elif test "${with_c4}" = mpi ; then
       test_flags="--${test_exe:=binary} --mpi"
   elif test "${with_c4}" = shmem ; then
       test_flags="--${test_exe:=binary} --shmem"
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

   dnl substitutions into Makefile.in
   AC_SUBST(NMLGEN)dnl
   AC_SUBST(LD)dnl
   AC_SUBST(AR)dnl
   AC_SUBST(ARFLAGS)dnl
   AC_SUBST(ARLIBS)dnl
   AC_SUBST(LIBDEPENDS)dnl

   # files to install

   : ${installfiles:='${install_executable} ${install_lib} ${install_headers}'}
   AC_SUBST(installfiles)dnl
   AC_SUBST(install_executable)dnl
   AC_SUBST(install_lib)dnl
   AC_SUBST(install_headers)dnl
   AC_SUBST(installdirs)dnl

   AC_SUBST(alltarget)dnl
   AC_SUBST(libsuffix)dnl
   AC_SUBST(dirstoclean)dnl
   AC_SUBST(package)dnl
   AC_SUBST(DRACO_DEPENDS)dnl
   AC_SUBST(DRACO_LIBS)dnl
   AC_SUBST(VENDOR_DEPENDS)dnl
   AC_SUBST(VENDOR_LIBS)dnl

   AC_SUBST(testcppflags)dnl
   AC_SUBST(PKG_DEPENDS)dnl
   AC_SUBST(PKG_LIBS)dnl
   AC_SUBST(DRACO_TEST_DEPENDS)dnl
   AC_SUBST(DRACO_TEST_LIBS)dnl
   AC_SUBST(VENDOR_TEST_DEPENDS)dnl
   AC_SUBST(VENDOR_TEST_LIBS)dnl
   AC_SUBST(test_alltarget)dnl
   AC_SUBST(test_flags)dnl
   AC_SUBST(test_scalar)dnl
   AC_SUBST(test_nprocs)dnl
   AC_SUBST(test_output_files)dnl

   AC_SUBST(configure_command)dnl

   CXXFLAGS="${CXXFLAGS} \${STRICTFLAG}"
   AC_SUBST(STRICTFLAG)

   dnl end of AC_DRACO_ENV
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_dracoenv.m4
dnl-------------------------------------------------------------------------dnl

