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
   dnl MAKE INSTALLS READ-ONLY
   dnl
       
   # no need to run AC_SUBST because it is done automagically by autoconf
   INSTALL_DATA="\${INSTALL} -m 444"

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
   dnl first find the system
   dnl

   AC_CANONICAL_SYSTEM

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

   dnl set up platform-dependent stuff in the 
   dnl SYSTEM-SPECIFIC SETUP SECTION
   
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
   
   dnl check for and assign the path to python
   AC_PATH_PROG(PYTHON_PATH, python, null)
   if test "${PYTHON_PATH}" = null ; then
       AC_MSG_ERROR("No valid Python found!")
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

   dnl
   dnl COMPILER SETUPS
   dnl

   dnl first find the system
   AC_CANONICAL_SYSTEM

   dnl set up a default compiler
   case $host in
   mips-sgi-irix6.*)   
       if test -z "${with_cxx}" ; then
	   with_cxx='cc'
       fi
   ;;
   alpha-dec-osf*)
       if test -z "${with_cxx}" ; then
	   with_cxx='compaq'
       fi
   ;;
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

   elif test "${with_cxx}" = guide ; then
       AC_CHECK_PROG(CXX, guidec++, guidec++)
       AC_CHECK_PROG(CC, guidec, guidec)

       if test "${CXX}" = guidec++ && test "${CC}" = guidec ; then 
	   AC_DRACO_GUIDE
       else
	   AC_MSG_ERROR("Did not find Guide compiler!")
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
   fi

   # do draco standard headers

   AC_MSG_CHECKING("for draco standard headers")
   if test "${enable_draco_stdhdrs:=no}" != no ; then
      CPPFLAGS="${CPPFLAGS} "'-I${includedir}/stdheaders'
      AC_MSG_RESULT("CPPFLAGS modified")
   else
      AC_MSG_RESULT("no") 
   fi

   # if with_f90 defined test with_f90 for compiler, and call setup
   # if with_f90 set to yes or not set 
   # attempt to guess compiler based on target

   if test "${with_f90+set}" = "set"   
   then
       AC_F90_ENV
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

   dnl determine and define data types for word sizes

   # eight byte integer types
   if test -n "${def_eight_byte_int_type}" ; then
       AC_DETERMINE_INT(8)
       AC_DEFINE_UNQUOTED(EIGHT_BYTE_INT_TYPE, ${INTEGER_SIZE_TYPE})
       if test "${INTEGER_SIZE_TYPE}" = 'long long' ; then
	   long_long_used='true'
       fi
   fi

   # four byte integer types
   if test -n "${def_four_byte_int_type}" ; then
       AC_DETERMINE_INT(4)
       AC_DEFINE_UNQUOTED(FOUR_BYTE_INT_TYPE, ${INTEGER_SIZE_TYPE})
       if test "${INTEGER_SIZE_TYPE}" = 'long long' ; then
	   long_long_used='true'
       fi
   fi

   # eight byte float types
   if test -n "${def_eight_byte_float_type}" ; then
       AC_DETERMINE_FLOAT(8)
       AC_DEFINE_UNQUOTED(EIGHT_BYTE_INT_TYPE, ${FLOAT_SIZE_TYPE})
   fi

   # four byte float types
   if test -n "${def_four_byte_float_type}" ; then
       AC_DETERMINE_FLOAT(4)
       AC_DEFINE_UNQUOTED(FOUR_BYTE_INT_TYPE, ${FLOAT_SIZE_TYPE})
   fi


   dnl
   dnl SYSTEM-SPECIFIC SETUP
   dnl

   # systems setup
   case $host in
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
       # setup communication packages
       #
   
       # setup for mpi support, on linux vendor and mich are one
       # and the same because there is no vendor for mpi on linux
        
       if test "${with_mpi}" = vendor ; then
	   with_mpi='mpich'
       fi

       if test "${with_mpi}" = mpich ; then
	   
	   # define mpich libraries for v1.2 of mpich
	   if test -n "${MPI_LIB}" ; then
	       linux_mpi_libs="-L${MPI_LIB} -lmpich"
	   elif test -z "${MPI_LIB}" ; then
	       linux_mpi_libs="-lmpich"
	   fi

	   # define the linux mpi libs
	   AC_VENDORLIB_SETUP(vendor_mpi, ${linux_mpi_libs})

       fi

       # shmem (not available on suns)
       if test "${enable_shmem}" = yes ; then
	   AC_MSG_ERROR("We do not support shmem on linux!")
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

	   # if an lapack location was defined use it
	   if test -n "${LAPACK_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_lapack, -L${LAPACK_LIB} -llapack -lblas)
	   elif test -z "${LAPACK_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_lapack, -llapack -lblas)
	   fi

       fi 

       # 
       # end of lapack setup
       # 

       #
       # setup eospac
       #
       
       AC_MSG_CHECKING("for extra eospac library requirements.")
       if test -n "${vendor_eospac}"; then
	   extra_eospac_libs="-L/usr/local/lf9560/lib -lfj9i6 -lfj9e6 -lfj9f6 -lfst -lfccx86"
           LIBS="${LIBS} ${extra_eospac_libs}"
           AC_MSG_RESULT("${extra_eospac_libs}")
       else
           AC_MSG_RESULT("none.")
       fi

       #
       # end of eospac
       #

       # add libg2c to LIBS if lapack, gandolf, or pcg is used
       AC_MSG_CHECKING("libg2c requirements")
       if test -n "${vendor_lapack}" || test -n "${vendor_pcg}" ||
	  test -n "${vendor_gandolf}"; then
	   
	   # if KCC compiler than add backend
	   if test "${CXX}" = KCC ; then
	       LIBS="${LIBS} --backend -lg2c"
	       AC_MSG_RESULT("--backend -lg2c added to LIBS")
	   elif test "${CXX}" = g++ ; then
	       LIBS="${LIBS} -lg2c"
	       AC_MSG_RESULT("-lg2c added to LIBS")
	   fi

       else
	   AC_MSG_RESULT("not needed")
       fi

       # set rpath when building shared library executables
       if test "${enable_shared}" = yes; then

	   # turn off ranlib
	   RANLIB=':'

	   # the g++ rpath needs Xlinker in front of it
	   if test "${CXX}" = g++; then
	       RPATHA="-Xlinker -rpath \${curdir}"
	       RPATHB="-Xlinker -rpath \${curdir}/.."
	       RPATHC="-Xlinker -rpath \${libdir}"
	       LDFLAGS="${RPATHA} ${RPATHB} ${RPATHC} ${LDFLAGS}"
	   else
	       LDFLAGS="-rpath \${curdir}:\${curdir}/..:\${libdir} ${LDFLAGS}"
	   fi

       fi
   ;;
   mips-sgi-irix6.*)
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
       elif test "${with_cxx}" = cc ; then
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
	   if test -n "${MPI_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_mpi, -L${MPI_LIB} -lmpi)
	   elif test -z "${MPI_LIB}" ; then
	       if test "${enable_32_bit}" = no ; then
		   AC_VENDORLIB_SETUP(vendor_mpi, -L${MPT_SGI}/usr/lib64 -lmpi)
	       else
		   AC_VENDORLIB_SETUP(vendor_mpi, -L${MPT_SGI}/usr/lib32 -lmpi)
	       fi
	   fi
       elif test "${with_mpi}" = mpich ; then
	   AC_MSG_ERROR("We do not support mpich on the SGI yet!")
       fi

       # setup for shmem support
       if test "${enable_shmem}" = yes ; then
	   if test -n "${SHMEM_LIB}" ; then
	       VENDOR_LIBS="${VENDOR_LIBS} -L${SHMEM_LIB} -lsma -lpthread"
	   elif test -z "${SHMEM_LIB}" ; then
	       VENDOR_LIBS="${VENDOR_LIBS} -lsma -lpthread"
	   fi
       fi

       # MPT (Message Passing Toolkit) for SGI vendor
       # implementation of MPI and SHMEM
       if test -z "${MPI_INC}" &&  test "${with_mpi}" = vendor ; then
	   MPI_INC="${MPT_SGI}/usr/include/"	   
	   MPI_H="\"${MPI_INC}mpi.h\""
	   AC_DEFINE_UNQUOTED(MPI_H, ${MPI_H})
       fi

       if test -z "${SHMEM_INC}" && test "${enable_shmem}" = yes ; then
	   SHMEM_INC="${MPT_SGI}/usr/include/mpp/"
	   SHMEM_H="\"${SHMEM_INC}shmem.h\""
	   AC_DEFINE_UNQUOTED(SHMEM_H, ${SHMEM_H})
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
	   if test "${with_mpi}" = vendor || test -n "${long_long_used}" ; then 
		   AC_MSG_WARN("KCC strict option set to allow long long type")
		   STRICTFLAG="${STRICTFLAG} --diag_suppress 450"
	   fi

       fi

       #
       # setup lapack
       #

       if test "${with_lapack}" = vendor ; then

	   # if an lapack location was defined use it
	   if test -n "${LAPACK_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_lapack, -L${LAPACK_LIB} -lcomplib.sgimath)
	   elif test -z "${LAPACK_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_lapack, -lcomplib.sgimath)
	   fi

       fi

       #
       # end of lapack setup
       #

       #
       # gandolf requires -lfortran on the link line.
       #

       AC_MSG_CHECKING("libfortran requirements")
       if test -n "${vendor_gandolf}" ; then
          LIBS="${LIBS} -lfortran"
          AC_MSG_RESULT("-lfortran added to LIBS")
       fi

       #
       # end of gandolf/libfortran setup
       #

       # set rpath when building shared library executables
       if test "${enable_shared}" = yes; then

	   # the g++ rpath needs Xlinker in front of it
	   if test "${CXX}" = g++; then
	       RPATHA="-Xlinker -rpath \${curdir}"
	       RPATHB="-Xlinker -rpath \${curdir}/.."
	       RPATHC="-Xlinker -rpath \${libdir}"
	       LDFLAGS="${RPATHA} ${RPATHB} ${RPATHC} ${LDFLAGS}"
	   else
	       LDFLAGS="-rpath \${curdir}:\${curdir}/..:\${libdir} ${LDFLAGS}"
	   fi

       fi
   ;;
   alpha-dec-osf*)
       # posix source defines, by default we set posix on 
       if test "${with_posix:=no}" = yes ; then
	   with_posix='199309L'
       fi

       if test "${with_posix}" != no ; then
	   AC_DEFINE_UNQUOTED(_POSIX_C_SOURCE, $with_posix)
	   AC_DEFINE(_POSIX_SOURCE)
       fi
   ;;
   sparc-sun-solaris2.*)
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
	   sun_mpi_libs='-lpmpi -lmpi -lsocket -lnsl'
   
	   if test -n "${MPI_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_mpi, -L${MPI_LIB} ${sun_mpi_libs})
	   elif test -z "${MPI_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_mpi, ${sun_mpi_libs})
	   fi

       fi

       # shmem (not available on suns)
       if test "${enable_shmem}" = yes ; then
	   AC_MSG_ERROR("We do not support shmem on suns!")
       fi

       #
       # end of communication package setup
       #

       #
       # setup lapack
       #

       if test "${with_lapack}" = vendor ; then

	   sun_libs='-llapack -lblas -lF77 -lM77 -lsunmath'

	   # if an lapack location was defined use it
	   if test -n "${LAPACK_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_lapack, -L${LAPACK_LIB} ${sun_libs})
	   elif test -z "${LAPACK_LIB}" ; then
	       AC_VENDORLIB_SETUP(vendor_lapack, ${sun_libs})
	   fi

       fi

       #
       # end of lapack setup
       #

       # set -R when building shared library executables
       if test "${enable_shared}" = yes; then
	   LDFLAGS="-R \${curdir}:\${curdir}/..:\${libdir} ${LDFLAGS}"
	   RANLIB=':'
       fi
   ;;
   *)
       # catchall for nothing
   ;;
   esac

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

   # compiler substitutions
   # automatic substitutions: CXXFLAGS, LDFLAGS, LIBS, CLFLAGS, CPPFLAGS
   AC_SUBST(LD)dnl
   AC_SUBST(AR)dnl
   AC_SUBST(ARFLAGS)dnl
   AC_SUBST(STRICTFLAG)dnl
   AC_SUBST(PARALLEL_FLAG)dnl

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
   AC_SUBST(VENDOR_LIBS)dnl
   AC_SUBST(ARLIBS)dnl

   # package testing libraries
   AC_SUBST(testcppflags)dnl
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

   # configure options
   AC_SUBST(configure_command)dnl

   # namelist package (nml) generation   
   AC_SUBST(NMLGEN)dnl

   dnl end of AC_DRACO_ENV
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_dracoenv.m4
dnl-------------------------------------------------------------------------dnl

