dnl-------------------------------------------------------------------------dnl
dnl ac_dracoarg.m4
dnl DRACO arguments macro that defines DRACO's non-vendor arguments
dnl
dnl $Id$
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl AC_DRACO_ARGS
dnl
dnl usage: configure.in
dnl defines non-vendor arguments for DRACO
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_ARGS, [dnl

   dnl
   dnl c4 toggle (scalar by default)
   dnl

   # define --with-c4
   AC_ARG_WITH(c4, 
      [  --with-c4[=scalar,mpi,shmem]   
		          turn on c4 (default scalar) ])

   # give with-c4 implied argument
   if test "${with_c4:=scalar}" = yes ; then
       with_c4='scalar'
   fi

   # now do the correct #defines
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

   # now set up the platform-independent directories
   AC_COMM_SET

   dnl
   dnl DBC toggle
   dnl

   dnl defines --with-dbc
   AC_ARG_WITH(dbc,
      [  --with-dbc[=level]      set Design-by-Contract])
	
   if test "${with_dbc}" = yes ; then
       with_dbc='7'
   elif test "${with_dbc}" = no ; then
       with_dbc='0'
   fi

   if test "${with_dbc:=default}" != default ; then
       AC_DEFINE_UNQUOTED(DBC, $with_dbc)
   fi
	
   dnl
   dnl SHARED versus ARCHIVE libraries
   dnl

   dnl defines --enable-shared
   AC_ARG_ENABLE(shared,
      [  --enable-shared         turn on shared libraries (.a default)])

   if test "${enable_shared:=no}" = yes ; then
       libsuffix='.so'
   else
       libsuffix='.a'
   fi

   dnl
   dnl STATIC VERSUS DYNAMIC LINKING
   dnl

   dnl defines --enable-static-ld
   AC_ARG_ENABLE(static-ld,
      [  --enable-static-ld      use (.a) libraries if possible])

   if test "${with_c4}" = shmem ; then
       enable_static_ld='yes'
   fi

   dnl
   dnl ANSI STRICT COMPLIANCE
   dnl

   dnl defines --enable-strict-ansi
   AC_ARG_ENABLE(strict-ansi,
      [  --disable-strict-ansi   turn off strict ansi compliance])

   dnl
   dnl COMPILER OPTIMZATION LEVEL
   dnl

   dnl defines --with-opt
   AC_ARG_WITH(opt,
      [  --with-opt[=0,1,2,3]    set optimization level (0 by default)])

   if test "${with_opt}" = yes ; then
       with_opt='0'
   fi

   dnl defines --enable-debug
   AC_ARG_ENABLE(debug,
      [  --enable-debug          turn on debug (-g) option])

   dnl
   dnl POSIX SOURCE
   dnl

   dnl defines --with-posix
   AC_ARG_WITH(posix,
      [  --with-posix[=num]      give posix source (199309L default)])

   if test "${with_posix}" = yes ; then
       AC_DEFINE(_POSIX_C_SOURCE, "199309L")
       AC_DEFINE(_POSIX_SOURCE)
   elif test "${with_posix:=199309L}" != no ; then
       AC_DEFINE_UNQUOTED(_POSIX_C_SOURCE, $with_posix)
       AC_DEFINE(_POSIX_SOURCE)
   fi

   dnl
   dnl CHOSE BIT COMPILATION ON SGI'S
   dnl

   dnl defines --enable-32-bit
   AC_ARG_ENABLE(32-bit,
      [  --enable-32-bit         do 32-bit compilation (SGI ONLY)])

   dnl
   dnl CHOSE MIPS INSTRUCTION SET ON SGI'S
   dnl

   dnl defines --with-mips
   AC_ARG_WITH(mips,
      [  --with-mips[=1,2,3,4]   set mips, mips4 by default (SGI ONLY)])

   if test "${with_mips}" = yes ; then
       with_mips='4'
   fi

   dnl end of AC_DRACO_ARGS
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_dracoarg.m4
dnl-------------------------------------------------------------------------dnl

