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
   LD='${CXX}'
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
dnl end of ac_compiler.m4
dnl-------------------------------------------------------------------------dnl

