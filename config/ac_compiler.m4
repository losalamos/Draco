dnl-------------------------------------------------------------------------dnl
dnl ac_compiler.m4
dnl sets up all of the compilers
dnl
dnl $Id$
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
   mips-sgi-irix6.*)   
       if test -z "${with_cxx}" ; then
	   with_cxx='sgi'
       fi
   ;;
   alpha*-dec-osf*)
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

   elif test "${with_cxx}" = icc ; then 
       AC_CHECK_PROG(CXX, icc, icc)

       if test "${CXX}" = icc ; then
	   CC='icc'
	   AC_DRACO_INTEL_ICC
       else
	   AC_MSG_ERROR("Did not find Intel icc compiler!")
       fi
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
dnl GUIDEC++ (KCC) COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_GUIDE, [dnl

   AC_MSG_CHECKING("configuration of ${CXX}/${CC} compilers")

   # KCC SPECIFIC FLAGS
   dirstoclean='ti_files'

   # LINKER AND LIBRARY (AR)
   LD='${CXX}'
   AR='${CXX}'
   ARFLAGS='--exceptions -o'
   ARLIBS='${DRACO_LIBS}'
   ARTESTLIBS='${PKG_LIBS} ${DRACO_TEST_LIBS} ${DRACO_LIBS}'

   # COMPILATION FLAGS

   # strict asci compliance
   if test "${enable_strict_ansi:=yes}" = yes ; then
       STRICTFLAG="--strict"
   fi

   # --one_per flag
   if test "${enable_one_per:=yes}" = yes ; then
       # yes there is an extra space before the flag
       ONEPERFLAG=" --one_per"
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
   # yes there is no space before the flag
   # also, guide turns off exceptions by default, here we turn them on

   CXXFLAGS="${CXXFLAGS}${ONEPERFLAG} --exceptions"

   # For version 3.3 of KCC the strict and thread_safe
   # cannot be used together (in general).

   if test "$with_c4" = shmem ; then
       CXXFLAGS="${CXXFLAGS} --thread_safe"
       STRICTFLAG=""
       LDFLAGS="${LDFLAGS} --thread_safe --static_libKCC"
   fi

   # add exceptions to ldflags
   LDFLAGS="--exceptions ${LDFLAGS}"

   AC_MSG_RESULT("Guide compiler flags set")

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
       STRICTFLAG="-ansi -Wnon-virtual-dtor"
   fi

   # optimization level
   # gcc allows -g with -O (like KCC)
    
   # defaults
   if test "${enable_debug:=yes}" = yes ; then
       gcc_debug_flag='-g'
   fi
   CXXFLAGS="${gcc_debug_flag} -O${with_opt:=0}"
   CFLAGS="${gcc_debug_flag} -O${with_opt:=0}"

   # add inlining if optimization is 01, 02 (it is on by default for
   # 03)
   if test "${with_opt}" = 1 || test "${with_opt}" = 2 ||
      test "${with_opt}" = 3 ; then
       CXXFLAGS="${CXXFLAGS} -finline-functions"
       CFLAGS="${CFLAGS} -finline-functions"
   fi

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

   dnl end of AC_DRACO_EGCS
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
       ARFLAGS='-shared -nocxxstd -expect_unresolved '*3td*' '
       ARFLAGS="${ARFLAGS} -expect_unresolved '*8_RWrwstd*' -o"
   else
       AR='ar'
       ARFLAGS='cr'
   fi

   # with the default template flags (-pt), the contents of the
   # cxx_repository do not seem to need adding when building
   # shared libraries; you do have to add them for archives
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
    
   # defaults
   if test "${enable_debug:=yes}" = yes ; then
       gcc_debug_flag='-g'
   fi
   CXXFLAGS="${CXXFLAGS} ${gcc_debug_flag} -O${with_opt:=0}"
   CFLAGS="${CFLAGS} ${gcc_debug_flag} -O${with_opt:=0}"

   # add ieee flag
   CXXFLAGS="${CXXFLAGS} -ieee"
   CFLAGS="${CFLAGS} -ieee"

   # turn off implicit inclusion
   CXXFLAGS="${CXXFLAGS} -noimplicit_include"

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
   # dirstoclean='ti_files'

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

   # --one_per flag
   #if test "${enable_one_per:=yes}" = yes ; then
   #    ONEPERFLAG="--one_per"
   #fi

   # optimization level
   if test "${enable_debug:=no}" = yes && \
      test "${with_opt:=0}" != 0 ; then
      CXXFLAGS="${CXXFLAGS} -g"
      CFLAGS="${CFLAGS} -g"
   fi
   CXXFLAGS="${CXXFLAGS} -O${with_opt:=0}"
   CFLAGS="${CFLAGS} -O${with_opt:=0}"

   # static linking option
   if test "${enable_static_ld}" = yes ; then
       LDFLAGS="${LDFLAGS} -static"
   fi

   # Parallel build flag
   
   #PARALLEL_FLAG="--parallel_build \${nj}"

   # final compiler additions
   #CXXFLAGS="${CXXFLAGS} ${ONEPERFLAG}"

   AC_MSG_RESULT("icc compiler flags set")
   
   dnl end of AC_DRACO_INTEL_ICC
])

dnl-------------------------------------------------------------------------dnl
dnl FORTRAN COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

builtin(include, ac_f90.m4)dnl AC_LANG_FORTRAN90 and AC_PROG_F90

dnl-------------------------------------------------------------------------dnl
dnl end of ac_compiler.m4
dnl-------------------------------------------------------------------------dnl
