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

           # KCC
           KCC)

               # setup linux strict if the compiler is KCC (also turn
               # off the warnings about long long being non-standard)
               AC_MSG_NOTICE([Linux KCC strict option set to allow long long type!])
               STRICTFLAG="--linux_strict -D__KAI_STRICT --diag_suppress 450"
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
       # add thread safety if we are using KCC on linux
       #
       if test "${CXX}" = KCC ; then
	   CFLAGS="--thread_safe ${CFLAGS}"
	   CXXFLAGS="--thread_safe ${CXXFLAGS}"
	   ARFLAGS="--thread_safe ${ARFLAGS}"
	   LDFLAGS="--thread_safe ${LDFLAGS}"
       fi

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

       #
       # setup eospac
       #
       
       AC_MSG_CHECKING("for extra eospac library requirements.")
       if test -n "${vendor_eospac}"; then
           lahey_lib_loc=`which lf95 | sed -e 's/bin\/lf95/lib/'`
	   extra_eospac_libs="-L${lahey_lib_loc} -lfj9i6 -lfj9e6 -lfj9f6 -lfst -lfccx86_6a"
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
           else
               AC_MSG_RESULT("not needed")
	   fi

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
       # finalize vendors
       #
       AC_VENDOR_FINALIZE

       AC_DBS_SETUP_RPATH(rpath)

       # add the intel math library for better performance when
       # compiling with intel
       if test "${CXX}" = icc; then
	   LIBS="$LIBS -limf"
       fi
])


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

       AC_DBS_SETUP_RPATH(rpath)

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

       if test "${with_lapack}" = vendor ; then
	   lapack_libs='-lcxmlp -lcxml'
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
       # gandolf, eospac, pcg, udm require -lfor on the link line.
       #

       AC_MSG_CHECKING("libfortran requirements")
       if test -n "${vendor_gandolf}" || test -n "${vendor_eospac}" ||
          test -n "${vendor_pcg}" || test -n "${vendor_udm}"; then
           LIBS="${LIBS} -lfor"
           AC_MSG_RESULT("-lfor added to LIBS")
       else
	   AC_MSG_RESULT("not needed")
       fi

       #
       # end of gandolf/libfortran setup
       #

       #
       # libpcg/libudm/libfmpi setup
       #

       AC_MSG_CHECKING("libfmpi requirements")
       if test -n "${vendor_pcg}" || test "${with_udm}" = mpi; then
           LIBS="${LIBS} -lfmpi"
           AC_MSG_RESULT("-lfmpi added to LIBS")
       else
	   AC_MSG_RESULT("not needed")
       fi

       #
       # end of libpcg setup
       #

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
       # finalize vendors
       #
       AC_VENDOR_FINALIZE

       AC_DBS_SETUP_RPATH(rpath)

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

       AC_DBS_SETUP_RPATH(R)
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

       AC_DBS_SETUP_RPATH(rpath)

]) dnl irix

dnl-------------------------------------------------------------------------dnl
dnl AC_DBS_DARWIN_ENVIRONMENT
dnl
dnl Configure draco build system Darwin-specific variables
dnl This function is called within AC_DBS_PLATFORM_ENVIRONMENT
dnl ***** NOT FULLY IMPLEMENTED *****
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_DARWIN_ENVIRONMENT], [dnl

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

       #
       # setup eospac
       #
       
       AC_MSG_CHECKING("for extra eospac library requirements.")
       if test -n "${vendor_eospac}"; then
           lahey_lib_loc=`which lf95 | sed -e 's/bin\/lf95/lib/'`
	   extra_eospac_libs="-L${lahey_lib_loc} -lfj9i6 -lfj9e6 -lfj9f6 -lfst -lfccx86_6a"
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
           else
               AC_MSG_RESULT("not needed")
	   fi

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
dnl $1 = rpath trigger.  One of "rpath" or "R"
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_SETUP_RPATH], [dnl

       rptrigger=$1
  
       if test "${enable_shared}" = yes ; then

	   # turn off ranlib
	   RANLIB=':'

	   # the g++/icc rpath needs Xlinker in front of it
	   if test "${CXX}" = g++ || test "${CXX}" = icc; then
	       RPATHA="-Xlinker -${rptrigger} \${curdir}"
	       RPATHB="-Xlinker -${rptrigger} \${curdir}/.."
	       RPATHC="-Xlinker -${rptrigger} \${libdir}"
	       RPATH="${RPATHA} ${RPATHB} ${RPATHC} ${RPATH}"
	   else
	       RPATH="-${rptrigger} \${curdir}:\${curdir}/..:\${libdir} ${RPATH}"
	   fi
       fi

       # add vendors to rpath
       for vendor_dir in ${VENDOR_LIB_DIRS}; 
       do
	   # if we are using gcc then add xlinker
	   if test "${CXX}" = g++ || test "${CXX}" = icc; then
	       RPATH="-Xlinker -${rptrigger} ${vendor_dir} ${RPATH}"

	   # else we just add the rpath
	   else
	       RPATH="-${rptrigger} ${vendor_dir} ${RPATH}"
	   fi
       done
]) dnl setup_rpath

dnl-------------------------------------------------------------------------dnl
dnl end of ac_platforms.m4
dnl-------------------------------------------------------------------------dnl
