dnl-------------------------------------------------------------------------dnl
dnl ac_platforms.m4
dnl
dnl Defines platform-specfic environments, including default vendor
dnl settings for the CCS-4/ASC computer platforms.
dnl
dnl Thomas M. Evans
dnl 2003/04/30 20:29:39
dnl-------------------------------------------------------------------------dnl
##---------------------------------------------------------------------------##
## $Id$
##---------------------------------------------------------------------------##

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

   # *****************************
   # MAC OS X/DARWIN POWERPC SETUP
   # *****************************
   powerpc-apple-darwin*)
       AC_DBS_DARWIN_PPC_ENVIRONMENT
   ;;      

   # *****************************
   # MAC OS X/DARWIN INTEL SETUP
   # *****************************
   *86-apple-darwin*)
       AC_DBS_DARWIN_INTEL_ENVIRONMENT
   ;;      

   # *******
   # NOTHING
   # *******
   *)
       AC_MSG_ERROR("Cannot figure out the platform or host!")
   ;;
   esac
])

dnl ------------------------------------------------------------------------ dnl
dnl AC_DBS_IFORT_ENVIRONMENT
dnl
dnl Some vendor setups require that the intel fortran compiler
dnl libraries be provided on the link line.  This m4 function adds the
dnl necessary libraries to LIBS.
dnl ------------------------------------------------------------------------ dnl
AC_DEFUN([AC_DBS_IFORT_ENVIRONMENT], [dnl

   # set the proper RPATH command depending on the C++ compiler
   case ${CXX} in 
       *g++ | *icpc | *ppu-g++)  rpath='-Xlinker -rpath ' ;;
       *pgCC)                    rpath='-R'               ;;
       *) AC_MSG_ERROR("Improper compiler set in LINUX.")
   esac

   AC_MSG_CHECKING("for extra ifort library requirements.")
   if test -n "${vendor_eospac}"    ||
      (test -n "${vendor_lapack}" && test "${with_lapack}" = "atlas") ||
      test -n "${vendor_scalapack}" ||
      test -n "${vendor_trilinos}"; then
      f90_lib_loc=`which ifort | sed -e 's/bin\/ifort/lib/'`
      extra_f90_libs="-L${f90_lib_loc} -lifcore -lgfortran"
      LIBS="${LIBS} ${extra_f90_libs}"
      AC_MSG_RESULT("${extra_f90_libs}")
   else
      AC_MSG_RESULT("none.")
   fi

   dnl Optimize flag   
   AC_MSG_CHECKING("for F90FLAGS")
   if test "${with_opt:=0}" != 0 ; then
      if test ${with_opt} -gt 2; then
         F90FLAGS="${F90FLAGS} -O3"
      else
         F90FLAGS="${F90FLAGS} -O${with_opt}"
      fi
   else 
      F90FLAGS="${F90FLAGS} -g"
   fi

   dnl C preprocessor flag
   F90FLAGS="${F90FLAGS}" 
   F90FREE="-cpp"
   AC_MSG_RESULT(${F90FLAGS})
   AC_MSG_RESULT(${F90FREE})

   dnl scalar or mpi ?
   AC_MSG_CHECKING("for F90MPIFLAGS")
   if test ${with_mpi:=no} = "no"; then
      F90FLAGS="${F90FLAGS} -DC4_SCALAR"
   else
       case ${with_mpi} in
       mpich)
         F90MPIFLAGS="-lfmpich"
         ;;
       lampi | LAMPI | LA-MPI)
         F90MPIFLAGS="-lmpi"
         ;;
       openmpi)
         F90MPIFLAGS ="-lmpi -lmpi_cxx -lmpi_f77" 
         ;;
       esac
   fi
   AC_MSG_RESULT(${F90MPIFLAGS})

   AC_MSG_CHECKING("for F90VENDOR_LIBS")
   F90VENDOR_LIBS="$F90VENDOR_LIBS ${F90MPIFLAGS} ${F90CXXFLAGS}"
   AC_MSG_RESULT("${F90VENDOR_LIBS}")
])

dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_LAHEY_ENVIRONMENT
dnl
dnl Some vendor setups require that the Lahey lib dir and compiler
dnl libraries be provided on the link line.  This m4 function adds the
dnl necessary libraries to LIBS.
dnl-------------------------------------------------------------------------dnl
AC_DEFUN([AC_DBS_LAHEY_ENVIRONMENT], [dnl

   if test `echo ${CXX} | sed -e 's/.*\///g'` != "g++"; then
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
   case `echo ${CXX} | sed -e 's/.*\///g'` in
       g++ | icpc | ppu-g++)
           rpath='-Xlinker -rpath '
           ;;
       pgCC)
           rpath='-R'
           ;;
       *)
           echo "Found CXX = $CXX.  Was looking for g++ | icpc | ppu-g++ | pgCC."
           AC_MSG_ERROR([Improper compiler set in LINUX (ac_dbs_pgf90_environment).])
   esac

   AC_MSG_CHECKING("for extra pgf90 library requirements.")
   if test -n "${vendor_eospac}"    ||
      (test -n "${vendor_lapack}" && test "${with_lapack}" = "atlas") ||
      test -n "${vendor_scalapack}" ||
      test -n "${vendor_trilinos}"; then
      f90_lib_loc=`which pgf90 | sed -e 's/bin\/pgf90/lib/'`
      # 64-bit pgf90 flags
      if test `uname -m` = x86_64 ; then
         extra_f90_libs="-L${f90_lib_loc}  -lpgf90rtl -lpgf90 -lpgf90_rpm1"
         extra_f90_libs="${extra_f90_libs}  -lpgf902 -lpgftnrtl -lpgc"
         extra_f90_rpaths="$rpath${f90_lib_loc}"
      else
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
   F90FLAGS="${F90FLAGS} -Mpreprocess -g77libs -pgcpplibs -pgf90libs"
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
   if test ${with_cxx} = "pgi"; then
      CXXLIBDIR=`which pgCC | sed -e 's/\/bin\/pgCC//'`
      CXXLIBDIR="${CXXLIBDIR}/lib"
      F90CXXFLAGS="-L${CXXLIBDIR} -lC -lstd"
   else
      CXXLIBDIR=${GCC_LIB_DIR}
dnl      F90CXXFLAGS="-L${CXXLIBDIR} -lstdc++"
      F90CXXFLAGS="-lrt"  
   fi   
   AC_MSG_RESULT(${F90CXXFLAGS})

   AC_MSG_CHECKING("for F90VENDOR_LIBS")
   F90VENDOR_LIBS="$F90VENDOR_LIBS ${F90MPIFLAGS} ${F90CXXFLAGS}"
   AC_MSG_RESULT("${F90VENDOR_LIBS}")
])

dnl ------------------------------------------------------------------------ dnl
dnl AC_DBS_GFORTRAN_ENVIRONMENT
dnl
dnl Some vendor setups require that the gfortran compiler libraries be provided 
dnl on the link line.  This m4 function adds the necessary libraries to LIBS.
dnl ------------------------------------------------------------------------ dnl
AC_DEFUN([AC_DBS_GFORTRAN_ENVIRONMENT], [dnl

   # set the proper RPATH command depending on the C++ compiler
   case `echo ${CXX} | sed -e 's/.*\///g'` in  
       g++ | ppu-g++)
           rpath='-Xlinker -rpath '
           ;;
       *)
           AC_MSG_ERROR("Improper compiler set in LINUX with gfortran.")
           ;;
   esac

   AC_MSG_CHECKING("for extra gfortran library requirements.")
   if test -n "${vendor_eospac}"    ||
      (test -n "${vendor_lapack}" && test "${with_lapack}" = "atlas") ||
      test -n "${vendor_scalapack}" ||
      test -n "${vendor_trilinos}"; then
      extra_f90_libs="-lgfortranbegin -lgfortran"
      LIBS="${LIBS} ${extra_f90_libs}"
      AC_MSG_RESULT("${extra_f90_libs}")
   else
      AC_MSG_RESULT("none.")
   fi

   dnl Optimize flag   
   AC_MSG_CHECKING("for F90FLAGS")
   if test "${with_opt:=0}" != 0 ; then
      if test ${with_opt} -gt 2; then
         F90FLAGS="${F90FLAGS} -O3"
      else
         F90FLAGS="${F90FLAGS} -O${with_opt}"
      fi
   else 
      F90FLAGS="${F90FLAGS} -g"
   fi

   dnl C preprocessor flag
   F90FLAGS="${F90FLAGS}" 
   F90FREE="-x f95-cpp-input"
   AC_MSG_RESULT(${F90FLAGS})
   AC_MSG_RESULT(${F90FREE})

   dnl scalar or mpi ?
   AC_MSG_CHECKING("for F90MPIFLAGS")
   if test ${with_mpi:=no} = "no"; then
      F90FLAGS="${F90FLAGS} -DC4_SCALAR"
   else
       case ${with_mpi} in
       mpich)
         F90MPIFLAGS="-lfmpich"
         ;;
       lampi | LAMPI | LA-MPI)
         F90MPIFLAGS="-lmpi"
         ;;
       openmpi)
         F90MPIFLAGS ="-lmpi -lmpi_cxx -lmpi_f77" 
         ;;
       esac
   fi
   AC_MSG_RESULT(${F90MPIFLAGS})

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
      test -n "${vendor_superludist}" ||
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

           case `echo ${CXX} | sed -e 's/.*\///g'` in

           # GNU g++
           g++) 
               AC_MSG_NOTICE([g++ -ansi option set to allow long long type!])
               STRICTFLAG="$STRICTFLAG -Wno-long-long"
           ;;

           # PGI
           pgCC)
               AC_MSG_NOTICE([pgCC suppressing diagnostic 450 to allow long long!])
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

       # setup F90 libs, rpath, etc. for apps when CXX is the
       # principal compiler
       if test "${with_f90:=no}" = no ; then
           case `echo ${CXX} | sed -e 's/.*\///g'` in
           pgCC)
               AC_CHECK_PROGS(F90, pgf90)
               ;;
           *)
               AC_CHECK_PROGS(F90, lf95 gfortran)
               ;;
           esac
           
           case ${F90} in
           *lf95)
               AC_MSG_CHECKING("if lahey found")
               AC_DBS_LAHEY_ENVIRONMENT
               ;;
           *pgf90)
               AC_DBS_PGF90_ENVIRONMENT
               ;;
           *gfortran)
               AC_DBS_GFORTRAN_ENVIRONMENT
               ;;
           *ifort)
               AC_DBS_IFORT_ENVIRONMENT
               ;;
           esac
       fi

       #
       # add librt to LIBS if udm is used
       #
       AC_MSG_CHECKING("librt requirements")
       if test -n "${vendor_udm}"; then

	   # Add rt for g++
           case $CXX in
           *g++)
	       LIBS="${LIBS} -lrt"
	       AC_MSG_RESULT("-lrt added to LIBS")
               ;;
           *)
               AC_MSG_RESULT("not needed")
               ;;
	   esac

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
               case $CXX in
               *g++)
                   CFLAGS="${CFLAGS} -fPIC"
                   CXXFLAGS="${CXXFLAGS} -fPIC"
                   AC_MSG_RESULT("-ldl added to LIBS -fPIC added to compile flags")
                   ;;
               *)
                   AC_MSG_RESULT("-ldl added to LIBS") ;;
               esac

           else  
               AC_MSG_RESULT("not needed")
           fi

       else
           AC_MSG_RESULT("not needed")
       fi

       #
       # PTHREAD FLAG: Add -pthread to CXXFLAGS if we are using either
       # Trilinos or STLPort
       #
       if test "${with_trilinos:-no}" != no || 
          test "${with_stlport:-no}" != no; then

          case ${CXX} in
            *g++) CXXFLAGS="${CXXFLAGS} -pthread" ;;
          esac
       fi

       #
       # Set up fpe_trap for this platform if gcc is on.
       #
       case ${CXX} in
         *g++)  AC_DEFINE(FPETRAP_LINUX_X86) ;;
       esac

       #
       # finalize vendors
       #
       AC_VENDOR_FINALIZE

       # handle rpaths
       if test "${with_f90:=no}" = no ; then
         case ${CXX} in
           *pgCC)         AC_DBS_SETUP_RPATH(-R, nospace) ;;
           *g++ | *icpc)  AC_DBS_SETUP_RPATH('-Xlinker -rpath', space) ;;
           *ppu-g++)      AC_DBS_SETUP_RPATH('-Xlinker -rpath', space) ;;
           *)             AC_MSG_ERROR("Unrecognized compiler on LINUX") ;;
         esac
       fi

       # add the intel math library for better performance when
       # compiling with intel
       case ${CXX} in
         *icpc) LIBS="$LIBS -limf" ;;
       esac

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
       # setup F90, libs, rpath, etc. for apps when CXX is the
       # principal compiler
       #
       if test "${with_f90:=no}" ; then
           AC_CHECK_PROGS(F90, f90)
           case ${F90} in
           f90)
               AC_DBS_COMPAQ_F90_ENVIRONMENT
               ;;
           esac
       fi

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
       if test "${with_cxx}" = ascipurple || test "${with_cxx}" = ibm; then
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
	   elif test "${with_cxx}" = ascipurple || 
                test "${with_cxx}" = ibm; then
	       CXXFLAGS="${CXXFLAGS} -q32"
	       CFLAGS="${CFLAGS} -q32"
	   fi

       elif test "${enable_64_bit:=no}" = yes ; then
	   
	   # switch on gcc or xlC compiler
	   if test "${with_cxx}" = gcc; then
	       CXXFLAGS="${CXXFLAGS} -maix64"
	       CFLAGS="${CFLAGS} -maix64"
	   elif test "${with_cxx}" = ascipurple || 
                test "${with_cxx}" = ibm; then
	       CXXFLAGS="${CXXFLAGS} -q64"
	       CFLAGS="${CFLAGS} -q64"
	   fi

       fi

       # set up the heap size
       if test "${with_cxx}" = ascipurple ; then
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

	       # on ascipurple the newmpxlC compiler script takes care
	       # of loading the mpi libraries; since it will fail
	       # if libraries are loaded and newmpxlC is used; throw
	       # an error if it occurs
	       if test "${with_cxx}" = ascipurple; then

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
	          test "${with_cxx}" = ascipurple ; then

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
dnl AC_DBS_DARWIN_COMMON_ENVIRONMENT
dnl
dnl Configure draco build system Darwin-specific variables, common to
dnl both PowerPC and Intel.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_DARWIN_COMMON_ENVIRONMENT], [dnl

       # dependency rules for IBM visual age compiler are complex
       if test "${with_cxx}" = ibm; then
	   DEPENDENCY_RULES='Makefile.dep.xlC.darwin'
       fi

       # print out cpu message
       AC_MSG_CHECKING("host platform cpu")
       AC_MSG_RESULT("${host_cpu}")

       AC_DBS_SETUP_POSIX

       #   
       # LONG LONG on Darwin
       #
       
       # always allow long long in strict ansi mode (if possible)
       
       if test -n "${STRICTFLAG}"; then

           case ${CXX} in

           # GNU g++
           *g++) 
               AC_MSG_NOTICE([g++ -ansi option set to allow long long type!])
               STRICTFLAG="$STRICTFLAG -Wno-long-long"
#               AC_MSG_NOTICE([g++ -ansi option set to allow long double type])
#               STRICTFLAG="$STRICTFLAG -Wno-long-double"
           ;;
  	   *ibm)	
	       AC_MSG_WARN("xlC set to allow long long")
	       STRICTFLAG="-qlanglvl=extended"
	       CFLAGS="${CFLAGS} -qlonglong"
	       CXXFLAGS="${CXXFLAGS} -qlonglong"
	   ;;

           *) # catchall
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
       AC_DBS_SETUP_COMM([${with_mpi:-openmpi}])

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

       # setup lf95 libs when CXX is the principle compiler
       if test "${with_f90:=no}" = no ; then
           AC_DBS_LAHEY_ENVIRONMENT
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
               case ${CXX} in
               *g++)
                   CFLAGS="${CFLAGS} -fPIC"
                   CXXFLAGS="${CXXFLAGS} -fPIC"
                   AC_MSG_RESULT("-ldl added to LIBS -fPIC added to compile flags")
                   ;;
               *)
                   AC_MSG_RESULT("-ldl added to LIBS") ;;
               esac

           else  
               AC_MSG_RESULT("not needed")
           fi

       else
           AC_MSG_RESULT("not needed")
       fi

])

dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_DARWIN_PPC_ENVIRONMENT
dnl
dnl Configure draco build system Darwin-PowerPC variables.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_DARWIN_PPC_ENVIRONMENT], [dnl

       AC_DBS_DARWIN_COMMON_ENVIRONMENT

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
dnl AC_DBS_DARWIN_INTEL_ENVIRONMENT
dnl
dnl Configure draco build system Darwin-Intel variables.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_DARWIN_INTEL_ENVIRONMENT], [dnl

       AC_DBS_DARWIN_COMMON_ENVIRONMENT

       #
       # Set up fpe_trap for this platform (not yet!!!)
       #
       AC_DEFINE(FPETRAP_DARWIN_INTEL)

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
        
       dnl echo "ac_dbs_setup_comm: 1        = $1"
       dnl echo "ac_dbs_setup_comm: with_mpi = $with_mpi"

       if test "${with_mpi}" = vendor ; then
	   with_mpi=$1
       fi

       # For CCS-2/4 users, we can also specify lampi or openmpi in place of
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
           if (expr " $mpi_version" : '.*LA-MPI' > /dev/null); then 
              AC_MSG_RESULT(${mpi_version})
           else
              AC_MSG_ERROR("Did not find LA-MPI version of mpirun.")
           fi
           ;;
       openmpi)
           with_mpi='OPENMPI'
           mpi_libs='-lmpi -lmpi_cxx -lmpi_f77'
           AC_MSG_CHECKING("mpirun -version")
           mpi_version=`mpirun -version 2>&1`
           if (expr " $mpi_version" : '.*Open MPI' > /dev/null); then 
              AC_MSG_RESULT(${mpi_version})
              # Do not include C++ bindings. See Draco artifact: artf7384
              CXXFLAGS="-DOMPI_SKIP_MPICXX ${CXXFLAGS}"
           else
              AC_MSG_ERROR("Did not find Open MPI version of mpirun.")
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
       for vendor_dir in ${VENDOR_LIB_DIRS}; do
           dnl Only append to RPATH if RPATH doesn't alreayd contain vendor_dir
           if (expr "${RPATH}" : ".*${vendor_dir}" > /dev/null); then
              continue
           fi
           dnl This is evil; should match against libdir, but libdir doesn't get
           dnl expanded until the actual make, so that's not an option here.
           if test "${vendor_dir}" = "${prefix}/lib"; then
              continue
           fi
           dnl Only append to RPATH if vendor has shared object libs.
           so_libs=`ls ${vendor_dir}/*.so 2>/dev/null`
           if test ! "${so_libs:-none}" = "none"; then
              if test "${dilem}" = "space"; then
	          RPATH="${rptrigger} ${vendor_dir} ${RPATH}"
              elif test "${dilem}" = "nospace"; then
	          RPATH="${rptrigger}${vendor_dir} ${RPATH}"
              elif test "${dilem}" = "colon"; then
	          RPATH="${rptrigger} ${vendor_dir} ${RPATH}"
              else
                  AC_MSG_ERROR("Cannot determine what rpath format to use!")
   	      fi
           fi 
       done

]) dnl setup_rpath

dnl-------------------------------------------------------------------------dnl
dnl end of ac_platforms.m4
dnl-------------------------------------------------------------------------dnl
