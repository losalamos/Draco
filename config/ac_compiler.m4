dnl-------------------------------------------------------------------------dnl
dnl ac_compiler.m4
dnl sets up all of the compilers
dnl
dnl $Id$
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl KCC COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_KCC, [dnl

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

   dnl end of AC_DRACO_KCC
])

dnl-------------------------------------------------------------------------dnl
dnl CC COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_CC, [dnl

   # LINKER AND LIBRARY (AR)
   LD='${CXX}'
   AR='${CXX}'
   ARLIBS='${DRACO_LIBS} ${VENDOR_LIBS}'

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

   dnl end of AC_DRACO_CC
])

dnl-------------------------------------------------------------------------dnl
dnl EGCS COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_EGCS, [dnl
   LD='${CXX}'
   dnl end of AC_DRACO_EGCS
])

dnl-------------------------------------------------------------------------dnl
dnl IBM XLF90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_XL_F90, [dnl

   # Check for working XL F90 compiler

   AC_MSG_CHECKING([for XLF90 compiler])
   AC_CHECK_PROG(F90, xlf90, xlf90, none)
   if test "${F90}" = xlf90
   then
       AC_MSG_RESULT([found])
   else
       AC_MSG_ERROR([not found])
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

   F90FLAGS="-qextchk -qhalt=s -qarch=pwr2 -bmaxstack:0x70000000 -bmaxdata:0x70000000 -qalias=noaryovrlp ${F90FREE}"

   if test "${enable_debug:=no}" = yes && test "${with_opt:=0}" != 0
   then
        trapflags="-qinitauto=FF"
        trapflags="${trapflags} -qflttrap=overflow:underflow:zerodivide:invalid:enable"
        trapflags="${trapflags} -qsigtrap"
        F90FLAGS="-g -d -C ${trapflags} -bloadmap:loadmap.dat ${F90FLAGS}"
   else
        F90FLAGS="-O${with_opt:=0} ${F90FLAGS}"
   fi

   dnl end of AC_DRACO_XL_F90
])

dnl-------------------------------------------------------------------------dnl
dnl FUJITSU F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_FUJITSU_F90, [dnl

   # Check for working Fujitsu F90 compiler

   AC_MSG_CHECKING([for Fujitsu f90 compiler])
   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90 && ${F90} -V 2>&1 | grep "Fujitsu"
   then
       AC_MSG_RESULT([found])
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

   # COMPILATION FLAGS

   F90FLAGS="-X9 -Am ${F90FREE}"

   if test "${enable_debug:=no}" = yes && test "${with_opt:=0}" != 0
   then
        F90FLAGS="-g -Haesu ${F90FLAGS}"
   else
        F90FLAGS="-O${with_opt:=0} ${F90FLAGS}"
   fi

   dnl end of AC_DRACO_FUJITSU_F90
])

dnl-------------------------------------------------------------------------dnl
dnl SUN WORKSHOP F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_WORKSHOP_F90, [dnl

   # Check for working WorkShop F90 compiler

   AC_MSG_CHECKING([for WorkShop f90 compiler])
   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90 && ${F90} -V 2>&1 | grep "WorkShop"
   then
       AC_MSG_RESULT([found])
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

   # Set COMPILATION FLAGS

   F90FLAGS="${F90FREE}"

   if test "${enable_debug:=no}" = yes && test "${with_opt:=0}" != 0
   then
        F90FLAGS="-g"
   else
        F90FLAGS="-O${with_opt:=0} ${F90FLAGS}"
   fi

   dnl end of AC_DRACO_WORKSHOP_F90
])

dnl-------------------------------------------------------------------------dnl
dnl CRAY_F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_CRAY_F90, [dnl

   # Check for working Cray F90 compiler

   AC_MSG_CHECKING([for Cray f90 compiler])
   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90
   then
       AC_MSG_RESULT([found])
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

   # COMPILATION FLAGS

   F90FLAGS="${F90FREE}"

   if test "${enable_debug:=no}" = yes && test "${with_opt:=0}" != 0
   then
        F90FLAGS="-g ${F90FLAGS}"
   else
        F90FLAGS="-O${with_opt:=0} ${F90FLAGS}"
   fi

   dnl end of AC_DRACO_CRAY_F90
])

dnl-------------------------------------------------------------------------dnl
dnl IRIX MIPS F90 COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_MIPS_F90, [dnl

   # Look for working MIPS compiler

   AC_MSG_CHECKING([for MIPS f90 compiler])
   AC_CHECK_PROG(F90, f90, f90, none)
   if test "${F90}" = f90 && ${F90} -version 2>&1 | grep "MIPS"
   then
       AC_MSG_RESULT([found])
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

   # COMPILATION FLAGS

   F90FLAGS="${F90FREE}"

   if test "${enable_debug:=no}" = yes && test "${with_opt:=0}" != 0
   then
        F90FLAGS="-g ${F90FLAGS}"
   else
        F90FLAGS="-O${with_opt:=0} ${F90FLAGS}"
   fi

   dnl end of AC_DRACO_MIPS_F90
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_compiler.m4
dnl-------------------------------------------------------------------------dnl
