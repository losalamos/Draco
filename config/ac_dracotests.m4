dnl-------------------------------------------------------------------------dnl
dnl ac_dracotests.m4 
dnl  
dnl Macros that run compiler/link tests after the draco environment
dnl has been configured.  These macros are called after AC_DRACO_ENV.
dnl
dnl Thomas M. Evans
dnl 2003/04/30 20:29:39
dnl-------------------------------------------------------------------------dnl
##---------------------------------------------------------------------------##
## $Id$
##---------------------------------------------------------------------------##

dnl-------------------------------------------------------------------------dnl
dnl AC_DETERMINE_WORD_TYPES
dnl
dnl defines the following variables that can be placed and used inside
dnl of a config.h file:
dnl
dnl   - SIZEOF_INT
dnl   - SIZEOF_LONG
dnl   - SIZEOF_LONG_LONG
dnl   - SIZEOF_FLOAT
dnl   - SIZEOF_DOUBLE
dnl   - SIZEOF_LONG_DOUBLE
dnl
dnl these are defined as macros, ie. #define SIZEOF_INT 4.  They can
dnl used to define types in a config.h file like so:
dnl
dnl     /* TYPES */
dnl     #undef SIZEOF_INT
dnl     #undef SIZEOF_LONG
dnl     #undef SIZEOF_LONG_LONG
dnl
dnl     /* SET THE TYPES */
dnl
dnl     /* Four byte int type */
dnl     #if SIZEOF_INT == 4
dnl     #define FOUR_BYTE_INT_TYPE int
dnl     #elif SIZEOF_LONG == 4
dnl     #define FOUR_BYTE_INT_TYPE long
dnl     #endif
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DETERMINE_WORD_TYPES], [dnl

   # the draco environment
   AC_REQUIRE([AC_DRACO_ENV])

   # integer types
   AC_CHECK_SIZEOF(int)
   AC_CHECK_SIZEOF(long)
   AC_CHECK_SIZEOF(long long)

   # float types
   AC_CHECK_SIZEOF(float)
   AC_CHECK_SIZEOF(double)
   AC_CHECK_SIZEOF(long double)
])

dnl-------------------------------------------------------------------------dnl
dnl AC_CHECK_EIGHT_BYTE_INT_TYPE
dnl
dnl This function should be called if the user plans on using an
dnl 8 byte integer type and wants to make sure that (1) an 8-byte type
dnl exists; and (2) any compiler strict flags are adjusted
dnl appropriately to use long long types
dnl
dnl If no 8-byte integer type is found a configure error is thrown
dnl
dnl This function should be called after AC_DETERMINE_WORD_TYPES
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_CHECK_EIGHT_BYTE_INT_TYPE], [dnl

   # make sure that the word types have been defined
   AC_REQUIRE([AC_DETERMINE_WORD_TYPES])

   # make sure that the host is defined
   AC_REQUIRE([AC_CANONICAL_HOST])

   # do we need to do mods if the compiler uses long long
   AC_MSG_CHECKING("8-byte integer type")

   # check on each int type to see if we need long long
   if test "${ac_cv_sizeof_int}" = 8; then
       AC_MSG_RESULT("int - no mods needed for $host/$CXX")
   elif test "${ac_cv_sizeof_long}" = 8; then
       AC_MSG_RESULT("long - no mods needed for $host/$CXX")
   elif test "${ac_cv_sizeof_long_long}" = 8; then
       
       # we have to do certain platform/compiler dependent mods to
       # use long long

       case $host in

       # *** SGI MODS
       mips-sgi-irix6.*)

	   # when we are using KCC we have to adjust the strict flag
	   # assuming it hasn't already been adjusted for the mpi
	   # vendor

	   if test "${CXX}" = KCC && 
	      test "${enable_strict_ansi}" = yes ; then
	   
	       # if integer type is long long and we aren't using the
	       # mpi vendor option---long long is already accounted
	       # for if the mpi vendor is on---then we need to adjust
	       # the strict flag
	       if test "${with_mpi}" != vendor; then 
		   AC_MSG_RESULT("KCC strict option set to allow long long type")
		   STRICTFLAG="${STRICTFLAG} --diag_suppress 450"
	       else
		   AC_MSG_RESULT("long long - no additional mods needed on $host/${CXX}")
	       fi

	   else

	       AC_MSG_RESULT("long long - no additional mods needed on $host/${CXX}")

	   fi
       ;;

       # *** IBM MODS
       *ibm-aix*)

	   # if we are using visual age then we may need to do some
	   # adjustment
	   if test "${with_cxx}" = ibm || 
	      test "${with_cxx}" = ascipurple ; then

	       # if the code package is serial we need to turn on long
	       # long or if mpi is on, but is not the vendor then
	       # we need long long
	       if test -z "${vendor_mpi}" || test "${with_mpi}" != vendor; then
	          
		   if test "${enable_strict_ansi}"; then
		       AC_MSG_RESULT("xlC set to allow long long")
		       STRICTFLAG="-qlanglvl=extended"
		       CFLAGS="${CFLAGS} -qlonglong"
		       CXXFLAGS="${CXXFLAGS} -qlonglong" 
		   else
		       AC_MSG_RESULT("long long - no additional mods needed on $host/${CXX}") 
		   fi
	   
	       else
		   AC_MSG_RESULT("long long - no additional mods needed on $host/${CXX}")
	       fi

	   else

	       AC_MSG_RESULT("long long - no additional mods needed on $host/${CXX}")
	   
	   fi
       ;;

       # *** OTHER SYSTEMS
       *)
	   AC_MSG_RESULT("long long - no additional mods needed on $host/${CXX}")
       ;;
       esac

   else
       AC_MSG_ERROR("no 8-byte int type available")
   fi
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_dracotests.m4
dnl-------------------------------------------------------------------------dnl
