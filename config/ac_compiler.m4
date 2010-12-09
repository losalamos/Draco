dnl-------------------------------------------------------------------------dnl
dnl ac_compiler.m4
dnl
dnl Sets up all of the C++ compilers.
dnl
dnl Thomas M. Evans
dnl 1999/03/05 18:16:55
dnl-------------------------------------------------------------------------dnl
##---------------------------------------------------------------------------##
## $Id$
##---------------------------------------------------------------------------##

dnl-------------------------------------------------------------------------dnl
dnl C++ COMPILER SETUP FUNCTION-->this is called within AC_DRACO_ENV;
dnl default is to use the C++ Compiler.  To change defaults,
dnl AC_WITH_F90 should be called in configure.in (before
dnl AC_DRACO_ENV)
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_CPP_ENV], [dnl

   # make sure that the host is defined
   AC_REQUIRE([AC_CANONICAL_HOST])

   dnl If not specified on the command line, set up a default compiler:
   dnl IRIX -> CC
   dnl COMPAQ -> CXX
   dnl IBM ASCI PURPLE -> newxlC (use --with-cxx=ibm for regular SP2)
   dnl EVERYTHING ELSE -> gcc
   if test -z "${with_cxx}" ; then
      case $host in
      mips-sgi-irix6.*) with_cxx='sgi' ;;
      alpha*-dec-osf*)  with_cxx='compaq' ;;
      *ibm-aix*)        with_cxx='ascipurple' ;;
      *)                with_cxx='gcc' ;;
      esac
   fi

   dnl determine which compiler we are using

   # do tests of --with-cxx, see if the compiler exists and then call
   # the proper setup function

   if test "${with_cxx}" = ppu-gcc || test "${with_cxx}" = ppu-g++ ; then 
       AC_CHECK_PROG(CXX, ppu-g++, ppu-g++)
       AC_CHECK_PROG(CC, ppu-gcc, ppu-gcc)

       if test `basename ${CXX}` = ppu-g++ && test `basename ${CC}` = ppu-gcc ; then
           AC_DRACO_GNU_PPU_GCC
       else
	   AC_MSG_ERROR("Did not find ppu-g++ compiler!")
       fi
   elif test "${with_cxx}" = sgi ; then
       AC_CHECK_PROG(CXX, CC, CC)
       AC_CHECK_PROG(CC, cc, cc)  

       if test `basename ${CXX}` = CC && test `basename ${CC}` = cc ; then
	   AC_DRACO_SGI_CC
       else 
	   AC_MSG_ERROR("Did not find SGI CC compiler!")
       fi

   elif test "${with_cxx}" = gcc ; then 
       AC_CHECK_PROG(CXX, g++, g++)
       AC_CHECK_PROG(CC, gcc, gcc)

       if test `basename ${CXX}` = g++ && test `basename ${CC}` = gcc ; then
	   AC_DRACO_GNU_GCC
       else
	   AC_MSG_ERROR("Did not find gnu c++ compiler!")
       fi

   elif test "${with_cxx}" = compaq ; then
       AC_CHECK_PROG(CXX, cxx, cxx)
       AC_CHECK_PROG(CC, cc, cc)

       if test `basename ${CXX}` = cxx && test `basename ${CC}` = cc ; then
	   AC_DRACO_COMPAQ_CXX
       else
	   AC_MSG_ERROR("Did not find Compaq cxx compiler!")
       fi

   elif test "${with_cxx}" = intel ; then 
       AC_CHECK_PROG(CXX, icpc, icpc)

       if test `basename ${CXX}` = icpc ; then
	   CC='icc'
	   AC_DRACO_INTEL_ICPC
       else
	   AC_MSG_ERROR("Did not find Intel icpc compiler!")
       fi

   elif test "${with_cxx}" = pgi ; then
       # only allow PGI on LINUX
       case $host in
       *-linux-gnu)
           AC_CHECK_PROG(CXX, pgCC, pgCC)
           
           if test `basename ${CXX}` == pgCC ; then 
               CC='pgcc'
               AC_DRACO_PGCC
           else
               AC_MSG_ERROR("Did not find PGI C++ compiler!")
           fi
       ;;
       *)
           AC_MSG_ERROR("PGI only available on LINUX.")
       ;;
       esac        

   elif test "${with_cxx}" = ibm ; then 
       AC_CHECK_PROG(CXX, xlC, xlC)
       AC_CHECK_PROG(CC, xlc, xlc)

       if test `basename ${CXX}` = xlC ; then
	   AC_DRACO_IBM_VISUAL_AGE
       else
	   AC_MSG_ERROR("Did not find IBM Visual Age xlC compiler!")
       fi

   elif test "${with_cxx}" = ascipurple ; then 
       # asci purple uses different executables depending upon the mpi
       # setup; so we check to see if mpi is on and set the executable
       # appropriately

       # mpi is on, use newmpxlC
       if test -n "${vendor_mpi}" && test "${with_mpi}" = vendor; then
	   AC_CHECK_PROG(CXX, newmpxlC, newmpxlC)
	   AC_CHECK_PROG(CC, newmpxlc, newmpxlc)

       # scalar build, use newxlC
       else
	   AC_CHECK_PROG(CXX, newxlC, newxlC)
	   AC_CHECK_PROG(CC, newxlc, newxlc)

       fi

       # check to make sure compiler is valid
       if test `basename ${CXX}` = newxlC || test `basename ${CXX}` = newmpxlC ; then
	   AC_DRACO_IBM_VISUAL_AGE
       else
	   AC_MSG_ERROR("Did not find ASCI Purple new(mp)xlC compiler!")
       fi

   else
       AC_MSG_ERROR("invalid compiler specification ${with_cxx}")

   fi

   # set the language to C++
   AC_LANG(C++)

])

dnl-------------------------------------------------------------------------dnl
dnl SGI CC COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DRACO_SGI_CC], [dnl

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

AC_DEFUN([AC_DRACO_GNU_GCC], [dnl

   # finding path of gcc compiler
   AC_PATH_PROG(GCC_BIN, g++, null)

   AC_MSG_CHECKING("Setting library path of GNU compiler")
   if test "${GCC_BIN}" = null ; then
       GCC_LIB_DIR='/usr/lib'
   else
       GCC_BIN=`dirname ${GCC_BIN}`
       GCC_HOME=`dirname ${GCC_BIN}`

       # Ensure that libraries exist at this location.  If we can't
       # find libstdc++.a at this location we leave GCC_LIB_DIR set to
       # null and issue a warning.

       dnl libstdc=`ls ${GCC_HOME}/lib/libstdc++.* | head -1`
       libstdc=`${CXX} -print-file-name=libstdc++.a`

       if test -n "${libstdc}" && test -r "${libstdc}"; then
         GCC_LIB_DIR="${GCC_HOME}/lib"
       fi
   fi
   AC_MSG_RESULT("${GCC_LIB_DIR}")

   if test -z ${GCC_LIB_DIR}; then
       AC_MSG_WARN("Could not determine location of gcc libraries. GCC_LIB_DIR is null")
   fi

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

   # Strict adhereance to the ISO C standard.
   dnl For C++ -ansi is equivalent to "-std=c++98"

   if test "${enable_strict_ansi:=yes}" = yes ; then
      STRICTFLAG="-ansi -Wnon-virtual-dtor -Wreturn-type -pedantic"
   fi

   # Verbosity of warnings
   # -Wall: This enables all the warnings about constructions that
   #    some users consider questionable, and that are easy to avoid
   #    (or modify to prevent the warning), even in conjunction with
   #    macros. 
   # -Wextra: This enables some extra warning flags that are not
   #    enabled by -Wall. 
   # -Weffc++: Warn about violations of the style guidelines from
   #    Scott Meyers' Effective C++ book. 
   # -Woverloaded-virtual: Warn when a function declaration hides
   #    virtual functions from a base class.
   # -Wcast-align: Warn whenever a pointer is cast such that the
   #    required alignment of the target is increased. 
   # -Wpointer-arith:  Warn about anything that depends on the "size
   #    of" a function type or of "void". 
   gcc_warn_flags="-Wall -Wextra -Weffc++ -Woverloaded-virtual "
   gcc_warn_flags="$gcc_warn_flags -Wcast-align -Wpointer-arith"

   # help for this variable activated in ac_dracoarg.m4:
   if test "${enable_all_warnings:=no}" = yes; then
      STRICTFLAG="${STRICTFLAG} ${gcc_warn_flags}"
   fi  

   # optimization level
   dnl Note: gcc allows -g with -O (like KCC)

   # set opt level in flags
   gcc_opt_flags="-O${with_opt:=0}"

   # set up compiler when optimized
   if test "${with_opt}" != 0; then

       # set up inlining when optimization is on
       gcc_opt_flags="-finline-functions ${gcc_opt_flags}"

       # turn off debug flag by default if not requested explicitly
       if test "${enable_debug:=no}" = yes ; then
	   gcc_opt_flags="-g ${gcc_opt_flags}"
       fi

   # set up compiler when not optimized
   else

       # Ref http://gcc.gnu.org/onlinedocs/libstdc++/manual/debug.html
       # default is to have debug flag on when opt=0
       if test "${enable_debug:=yes}" = yes ; then
           addflags="-g -fno-inline -fno-eliminate-unused-debug-types"
	   gcc_opt_flags="${addflags} ${gcc_opt_flags}"
           if test "${enable_glibcxx_debug:=no}" = yes; then
              addflags="-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC"
              gcc_opt_flags="${gcc_opt_flags} ${addflags}"
           fi
       fi

   fi

   # 64-bit gcc require -fPIC for enabled shared
   if test "${enable_shared}" = yes && test `uname -m` = x86_64 ; then
      CXXFLAGS="${CXXFLAGS} -fPIC"
      CFLAGS="${CFLAGS} -fPIC"
   fi
   
   # add opt flags
   CXXFLAGS="${gcc_opt_flags} ${CXXFLAGS}"
   CFLAGS="${gcc_opt_flags} ${CFLAGS}"

   # RPATH FLAGS

   # add -rpath for the compiler library (G++ as LD does not do this
   # automatically) if required.
   case $host in

   # Darwin doesn't need any special flags
   *-apple-darwin*)
   ;;

   # COMPAQ -> CXX
   alpha*-dec-osf*)
   ;;

   # EVERYTHING ELSE -> linux?
   *)
      if test -n "${GCC_LIB_DIR}"; then
           RPATH="${RPATH} -Xlinker -rpath ${GCC_LIB_DIR}"
      fi
   ;;
   esac

   # static linking option
   if test "${enable_static_ld}" = yes ; then
       LDFLAGS="${LDFLAGS} -Bstatic"
   fi

   AC_MSG_RESULT("GNU g++ compiler flags set")

   dnl end of AC_DRACO_GNU_GCC
])

dnl-------------------------------------------------------------------------dnl
dnl PPU GNU COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DRACO_GNU_PPU_GCC], [dnl

   # finding path of gcc compiler
   AC_PATH_PROG(PPU_GCC_BIN, ppu-g++, null)

   AC_MSG_CHECKING("Setting library path of ppu-g++ compiler")
   if test "${PPU_GCC_BIN}" = null ; then
       #PPU_GCC_LIB_DIR='/usr/lib'
       PPU_GCC_LIB_DIR='/opt/cell/toolchain/lib/gcc/ppu/4.1.1/'
   else
       PPU_GCC_BIN=`dirname ${PPU_GCC_BIN}`
       PPU_GCC_HOME=`dirname ${PPU_GCC_BIN}`

       # Ensure that libraries exist at this location.  If we can't
       # find libstdc++.a at this location we leave PPU_GCC_LIB_DIR set to
       # null and issue a warning.

       #if test -r ${PPU_GCC_HOME}/lib/libstdc++.a; then
       #  PPU_GCC_LIB_DIR="${PPU_GCC_HOME}/lib"
       #fi
       PPU_GCC_LIB_DIR='/opt/cell/toolchain/lib/gcc/ppu/4.1.1/'
   fi
   AC_MSG_RESULT("${PPU_GCC_LIB_DIR}")

   if test -z ${PPU_GCC_LIB_DIR}; then
       AC_MSG_WARN("Could not determine location of ppu_gcc libraries. PPU_GCC_LIB_DIR is null")
   fi

   # do compiler configuration
   AC_MSG_CHECKING("configuration of ${CXX}/${CC} compilers")

   # LINKER AND LIBRARY (AR)
   LD='${CXX}'

   # if shared then ar is ppu_gcc
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
       STRICTFLAG="-ansi -Wnon-virtual-dtor -Wreturn-type -pedantic"
   fi

   # optimization level
   # ppu_gcc allows -g with -O (like KCC)

   # set opt level in flags
   ppu_gcc_opt_flags="-O${with_opt:=0}"

   # set up compiler when optimized
   if test "${with_opt}" != 0; then

       # set up inlining when optimization is on
       ppu_gcc_opt_flags="-finline-functions ${ppu_gcc_opt_flags}"

       # turn off debug flag by default if not requested explicitly
       if test "${enable_debug:=no}" = yes ; then
	   ppu_gcc_opt_flags="-g ${ppu_gcc_opt_flags}"
       fi

   # set up compiler when not optimized
   else

       # default is to have debug flag on when opt=0
       if test "${enable_debug:=yes}" = yes ; then
	   ppu_gcc_opt_flags="-g ${ppu_gcc_opt_flags}"
       fi

   fi

   # 64-bit ppu_gcc require -fPIC for enabled shared
   if test "${enable_shared}" = yes && test `uname -m` = x86_64 ; then
      CXXFLAGS="${CXXFLAGS} -fPIC"
   fi
   
   # add opt flags
   CXXFLAGS="${ppu_gcc_opt_flags} ${CXXFLAGS}"
   CFLAGS="${ppu_gcc_opt_flags} ${CFLAGS}"

   # RPATH FLAGS

   # add -rpath for the compiler library (G++ as LD does not do this
   # automatically) if required.
   case $host in

   # Darwin doesn't need any special flags
   *-apple-darwin*)
   ;;

   # COMPAQ -> CXX
   alpha*-dec-osf*)
   ;;

   # EVERYTHING ELSE -> linux?
   *)
      if test -n "${PPU_GCC_LIB_DIR}"; then
           RPATH="${RPATH} -Xlinker -rpath ${PPU_GCC_LIB_DIR}"
      fi
   ;;
   esac

   # static linking option
   if test "${enable_static_ld}" = yes ; then
       LDFLAGS="${LDFLAGS} -Bstatic"
   fi

   AC_MSG_RESULT("ppu_g++ compiler flags set")

   dnl end of AC_DRACO_GNU_PPU_GCC
])

dnl-------------------------------------------------------------------------dnl
dnl PGI COMPILER SETUP
dnl 
dnl Note that this implementation of PGI uses options that are only
dnl valid for LINUX
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DRACO_PGCC], [dnl

   # do compiler configuration
   AC_MSG_CHECKING("configuration of ${CXX}/${CC} compilers")

   # LINKER AND LIBRARY (AR)
   LD='${CXX}'

   # if shared then ar is pgCC
   if test "${enable_shared}" = yes ; then
       AR="${CXX}"
       ARFLAGS='-fPIC -shared -o'
       LDFLAGS='-fPIC'

       # must use position-independent code
       CXXFLAGS="${CXXFLAGS} -fPIC"
       CFLAGS="${CFLAGS} -fPIC"
   else
       AR='ar'
       ARFLAGS='cr'
   fi

   ARLIBS=''
   ARTESTLIBS=''

   # COMPILATION FLAGS

   # strict asci compliance
   if test "${enable_strict_ansi:=yes}" = yes ; then
       STRICTFLAG="-Xa -A --no_using_std"

       # pgCC 9 and 10 have problem with our redhat systems.  
       # http://bit.ly/az9QIa
       # The suggested work around is to add -DNO_PGI_OFFSET to the
       # compile flags
       STRICTFLAG="${STRICTFLAG} -DNO_PGI_OFFSET"

       # suppress long long errors in the platform-dependent options
       # section 

       # suppress missing return statement warning (we get this in
       # nearly every STL inclusion through PGICC)
       STRICTFLAG="--diag_suppress 940 ${STRICTFLAG}"

       # suppress "unrecognized preprocessing directive" when pgCC
       # encounters #warning; because #warning is non-standard, emitting
       # a warning via #warning won't work anyway in strict mode
       STRICTFLAG="--diag_suppress 11 ${STRICTFLAG}"
   fi

   # optimization level
   # pgCC allows -g with -O

   # set opt level in flags

   # Consider adding: 
   # -Mipa=fast      invoke interprocedural analysis.
   # -Minline=levels:10
   # --no_exceptions
   # -Mpfi
   # -Mpfo
   # -Msafeptr
   # -O[0-4]
   pgcc_opt_flags="-O${with_opt:=0}"

   # set up compiler when optimized
   if test "${with_opt}" != 0; then

       # set up inlining when optimization is on
       pgcc_opt_flags="${pgcc_opt_flags}"

       # turn off debug flag by default if not requested explicitly
       if test "${enable_debug:=no}" = yes ; then
	   pgcc_opt_flags="-g ${pgcc_opt_flags}"
       fi

   # set up compiler when not optimized
   else

       # default is to have debug flag on when opt=0
       if test "${enable_debug:=yes}" = yes ; then
	   pgcc_opt_flags="-g ${pgcc_opt_flags}"
       fi
 
       # -c    array bounds checking
   fi

   # add opt flags
   CXXFLAGS="${pgcc_opt_flags} ${CXXFLAGS}"
   CFLAGS="${pgcc_opt_flags} ${CFLAGS}"
   
   # add ieee flag
   CXXFLAGS="${CXXFLAGS} -Kieee"
   CFLAGS="${CFLAGS} -Kieee"

   # instantiate only functions that are used in the compilation
   CXXFLAGS="${CXXFLAGS} --no_implicit_include"

   # set unnormalized values to zero
   CXXFLAGS="${CXXFLAGS} -Mdaz"
   CFLAGS="${CFLAGS} -Mdaz"

   AC_MSG_RESULT("PGI pgCC compiler flags set")

   dnl end of AC_DRACO_PGCC
])

dnl-------------------------------------------------------------------------dnl
dnl COMPAQ CXX COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DRACO_COMPAQ_CXX], [dnl

   dnl 6-FEB-02 NEED TO ADD MODS !!!!!!

   AC_MSG_CHECKING("configuration of ${CXX}/${CC} compilers")

   # CXX SPECIFIC FLAGS
   dirstoclean='cxx_repository'

   # LINKER AND LIBRARY (AR)
   LD='${CXX}'

   # if shared then ar is cxx
   if test "${enable_shared}" = yes ; then
       AR='${CXX}'
       ARFLAGS="-shared -nocxxstd"
       ARFLAGS="${ARFLAGS} -o"
   else
       AR='ar'
       ARFLAGS='cr'
   fi

   # the contents of the cxx_repository do not seem to need adding 
   # when building shared libraries; you do have to add them for
   # archives 
   if test "${enable_shared}" != yes ; then
       ARLIBS='$(wildcard cxx_repository/*)'
       ARTESTLIBS='$(wildcard cxx_repository/*)'
   fi

   # COMPILATION FLAGS

   # strict asci compliance
   if test "${enable_strict_ansi:=yes}" = yes ; then
       STRICTFLAG="-std strict_ansi"
       CXX="${CXX} -model ansi"
   fi

   # make sure we always use the standard IO stream
   CPPFLAGS="${CPPFLAGS} -D__USE_STD_IOSTREAM" 

   # optimization level

   # if optimization is on turn off debug flag unless asked for
   if test "${with_opt:=0}" != 0 ; then

       # if debug is on then use -g1,2,3
       if test "${enable_debug:=no}" = yes ; then
	   cxx_opt_flag="-g${with_opt}"
       else
	   cxx_opt_flag="-O${with_opt}"
       fi

   # turn off optimizations
   else
   
       # we want -g unless not asked for
       if test "${enable_debug:=yes}" = yes ; then
	   cxx_opt_flag="-g -O0"
       else
	   cxx_opt_flag="-O0"
       fi

   fi

   # set up cxx flags
   CXXFLAGS="${CXXFLAGS} ${cxx_opt_flag}"
   CFLAGS="${CFLAGS} ${cxx_opt_flag}"

   # add ieee flag
   CXXFLAGS="${CXXFLAGS} -ieee"
   CFLAGS="${CFLAGS} -ieee"

   # turn off implicit inclusion
   CXXFLAGS="${CXXFLAGS} -noimplicit_include"

   # use the -pt template option for the compiler:
   # -pt Automatically instantiate templates into the repository with
   #  external linkage. Manually instantiated templates are placed in
   #  the output object with external linkage. This option is the default.
   CXXFLAGS="${CXXFLAGS} -pt"

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
dnl Intel icpc COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DRACO_INTEL_ICPC], [dnl

   AC_MSG_CHECKING("configuration of ${CXX}/${CC} compilers")

   # icpc SPECIFIC FLAGS

   # LINKER AND LIBRARY
   LD='${CXX}'

   # if shared then ar is icpc
   if test "${enable_shared}" = yes ; then
       AR="${CXX}"
       ARFLAGS='-fPIC -shared -o'
       LDFLAGS='-fPIC'

       # must use position-independent code
       CXXFLAGS="${CXXFLAGS} -fPIC"
       CFLAGS="${CFLAGS} -fPIC"
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

   # set up compiler when optimized (enable inline keyword but not
   # compiler-choice inlining)
   if test "${with_opt:=0}" != 0 ; then

       # turn off debug by default
       if test "${enable_debug:=no}" = yes ; then
	   icpc_opt_flags="-g -O${with_opt} -inline-level=1 -ip"
           # icpc 10.0.023 apparently can't link when opt>0
           LDFLAGS="${LDFLAGS} -O0"
       else
	   icpc_opt_flags="-O${with_opt} -inline-level=1"
           # icpc 10.0.023 apparently can't link when opt>0
           LDFLAGS="${LDFLAGS} -O0"
       fi

   #set up compiler when not optimized (turn off inlining with -Ob0)
   else

       # turn on debug by default
       if test "${enable_debug:=yes}" = yes ; then
	   icpc_opt_flags="-g -O0 -inline-level=0"
           LDFLAGS="${LDFLAGS} -O0"
       else
	   icpc_opt_flags="-O0 -inline-level=0"
           LDFLAGS="${LDFLAGS} -O0"
       fi

   fi
   
   # set the cxx and c flags
   CXXFLAGS="${CXXFLAGS} ${icpc_opt_flags}"
   CFLAGS="${CFLAGS} ${icpc_opt_flags}"

   # static linking option
   if test "${enable_static_ld}" = yes ; then
       LDFLAGS="${LDFLAGS} -static"
   fi

   AC_MSG_RESULT("icpc compiler flags set")
   
   dnl end of AC_DRACO_INTEL_ICPC
])

dnl-------------------------------------------------------------------------dnl
dnl IBM VISUAL AGE COMPILER SETUP
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DRACO_IBM_VISUAL_AGE], [dnl

   AC_MSG_CHECKING("configuration of ${CXX}/${CC} compilers")

   # xlC SPECIFIC FLAGS

   # LINKER AND LIBRARY
   LD='${CXX}'

   # if shared then ar is xlC
   if test "${enable_shared}" = yes ; then
       AR="${CXX}"
       ARFLAGS='-brtl -Wl,-bh:5 -G -o'

       # when AR=newmpxlC we need to add /lib/crt0_64.o to 
       # avoid p_argcx and p_argvx link error when building libs
       if test "${AR}" = newmpxlC ; then
	   ARLIBS='/lib/crt0_64.o'
	   ARTESTLIBS='/lib/crt0_64.o'
       fi

       ARLIBS="${ARLIBS} \${DRACO_LIBS} \${VENDOR_LIBS}"
       ARTESTLIBS="${ARTESTLIBS} \${PKG_LIBS} \${DRACO_TEST_LIBS}"
       ARTESTLIBS="${ARTESTLIBS} \${DRACO_LIBS}\${VENDOR_TEST_LIBS}"
       ARTESTLIBS="${ARTESTLIBS} \${VENDOR_LIBS}"
   else
       AR='ar'
       ARFLAGS='cr'

       ARLIBS=''
       ARTESTLIBS=''
   fi

   # COMPILATION FLAGS

   # strict asci compliance
   if test "${enable_strict_ansi:=yes}" = yes ; then
       STRICTFLAG="-qlanglvl=strict98"
   fi

   # the qinline option controls inlining, when -g is on no inlining
   # is done, with -O# inlining is on by default

   # set up compiler when optimized 
   if test "${with_opt:=0}" != 0; then

       # optflags
       xlC_opt_flags="-qarch=auto -qtune=auto -qcache=auto"

       # optimization level    
       if test "${with_opt}" = 1; then
	   # if asking for 1 just use opt in ibm   
	   xlC_opt_flags="${xlC_opt_flags} -qopt"
       else
	   # otherwise use number

	   # turn of aggressive semantic optimizations on all levels
	   # -O2 and above
	   xlC_opt_flags="${xlC_opt_flags} -qopt=${with_opt} -qstrict"
       fi

       # turn off debug by default
       if test "${enable_debug:=no}" = yes ; then
	   xlC_opt_flags="-g ${xlC_opt_flags}"
       fi

   #set up compiler when not optimized 
   else

       # optflags
       xlC_opt_flags="-qnoopt"

       # turn on debug by default
       if test "${enable_debug:=yes}" = yes ; then
	   xlC_opt_flags="-g ${xlC_opt_flags}"
       fi

   fi
   
   # set the CXX and CC flags

   # set the optimizations
   CXXFLAGS="${CXXFLAGS} ${xlC_opt_flags}"
   CFLAGS="${CFLAGS} ${xlC_opt_flags}"

   # set template stuff
   CXXFLAGS="${CXXFLAGS} -w -qnotempinc"

   # static linking option
   if test "${enable_static_ld:=no}" = yes ; then
       LDFLAGS="${LDFLAGS} -bstatic"

   # if we are building shared libraries we need to add
   # run-time-linking
   else
       LDFLAGS="${LDFLAGS} -brtl -Wl,-bh:5"

   fi

   AC_MSG_RESULT("${CXX} compiler flags set")
   
   dnl end of AC_DRACO_IBM_VISUAL_AGE
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_compiler.m4
dnl-------------------------------------------------------------------------dnl
