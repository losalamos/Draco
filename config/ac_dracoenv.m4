dnl-------------------------------------------------------------------------dnl
dnl ac_dracoenv.m4
dnl puts together the DRACO environments given the arguments from 
dnl vendors (ac_vendors.m4) and DRACO (ac_dracoarg.m4)
dnl
dnl Time-stamp: <99/02/24 10:34:01 rsqrd>
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
   dnl COMPILER SETUPS
   dnl

   dnl first find the host
   AC_CANONICAL_HOST
   
   dnl test for the presence of KCC
   AC_CHECK_PROG(CXX, KCC, KCC, NO_KCC)
   AC_CHECK_PROG(CC, KCC, KCC --c, NO_KCC)
   if test "$CXX" = NO_KCC ; then
       AC_PROG_CXX
   else
       # KCC SPECIFIC FLAGS
       dirstoclean='ti_files'

       # LINKER AND LIBRARY (AR)
       LD='${CXX}'
       AR='${CXX}'
       ARFLAGS='-o'
       ARLIBS='${DRACO_LIBS} ${VENDOR_LIBS}'

       # COMPILATION FLAGS

       # strict asci compliance
       if test "${enable_strict_ansi:=yes}" = yes ; then
	   STRICTFLAG="--strict"
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

       # final compiler additions
       CXXFLAGS="${CXXFLAGS} --one_per"

       # For version 3.3 of KCC the strict and thread_safe
       # cannot be used together (in general).

       if test "$with_c4" = shmem ; then
          CXXFLAGS="${CXXFLAGS} --thread_safe"
          STRICTFLAG=""
          LDFLAGS="${LDFLAGS} --thread_safe --static_libKCC"
       fi
   fi

   if test "$CC" = NO_KCC ; then
       AC_PROG_CC
   fi
	
   dnl check for ranlib
   AC_PROG_RANLIB

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
       CXXFLAGS="-mips${with_mips:=4} --backend -r10000 ${CXXFLAGS}"
       CFLAGS="-mips${with_mips:=4} -r10000 ${CFLAGS}"
       LDFLAGS="-mips${with_mips:=4} ${LDFLAGS}"

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
   ;;
   *)
   esac

   # add system specific libraries
   LIBS='-lm'

   dnl
   dnl DEJAGNU TEST SYSTEM
   dnl

   # Double dollar signs in MAKE commands will become single dollar signs.
   # Backslashes are needed to ensure that objdir is not
   # evaluated until unix.exp

   for tool in $test_alltarget; do
       if test "${with_c4:=scalar}" = scalar || \
	   test "$with_c4" = yes ; then
	   test_launch='\$$objdir/'"${tool}"
	   test_nprocs="1"
       elif test "$with_c4" = mpi ; then
	   test_launch='mpirun -np \$$NPROCS \$$objdir/'"${tool}"
       elif test "$with_c4" = shmem ; then
	   test_launch='\$$objdir/'"${tool}"' -npes \$$NPROCS'
       fi

       test_output_files="$test_output_files $tool.sum $tool.log"
       site_exp="$site_exp 'set ${tool}name \"${test_launch}\"'"
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

   : ${RUNTEST:=runtest}
   : ${RUNTESTTOOLFLAG:="--tool"}
   : ${SITE_EXP:=site.exp}
   : ${test_nprocs:="1 2 4 7"}

   AC_SUBST(testcppflags)dnl
   AC_SUBST(PKG_DEPENDS)dnl
   AC_SUBST(PKG_LIBS)dnl
   AC_SUBST(DRACO_TEST_DEPENDS)dnl
   AC_SUBST(DRACO_TEST_LIBS)dnl
   AC_SUBST(VENDOR_TEST_DEPENDS)dnl
   AC_SUBST(VENDOR_TEST_LIBS)dnl
   AC_SUBST(test_alltarget)dnl
   AC_SUBST(RUNTEST)dnl
   AC_SUBST(RUNTESTFLAGS)dnl
   AC_SUBST(RUNTESTTOOLFLAG)dnl
   AC_SUBST(FLAGS_TO_PASS)dnl
   AC_SUBST(test_launch)dnl
   AC_SUBST(test_nprocs)dnl
   AC_SUBST(SITE_EXP)dnl
   AC_SUBST(site_exp)dnl
   AC_SUBST(test_output_files)dnl

   AC_SUBST(configure_command)dnl

   CXXFLAGS="${CXXFLAGS} \${STRICTFLAG}"
   AC_SUBST(STRICTFLAG)

   dnl end of AC_DRACO_ENV
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_dracoenv.m4
dnl-------------------------------------------------------------------------dnl

