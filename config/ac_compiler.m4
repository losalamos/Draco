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

   CXXFLAGS="${CXXFLAGS}${ONEPERFLAG}"

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
dnl end of ac_compiler.m4
dnl-------------------------------------------------------------------------dnl
