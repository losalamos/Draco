dnl-------------------------------------------------------------------------dnl
dnl ac_draco_vendor.m4
dnl macros for setting up draco as a vendor
dnl
dnl $Id$
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_VENDOR_SETUP, [dnl

   dnl define --with-draco 
   AC_ARG_WITH(draco,
      [  --with-draco     give head location of draco includes and libs])

   AC_ARG_WITH(draco-inc,
      [  --with-draco-lib give location of draco includes])

   AC_ARG_WITH(draco-lib,
      [  --with-draco-inc give location of draco libs])

   # draco include and lib dirs
   DRACO_INC=''
   DRACO_LIB=''

   AC_MSG_CHECKING("checking draco directories") 

   # if we use the with-draco option to specify the location of draco,
   # then check it
   if test -n "${with_draco}" ; then 

       # check to see if we should set draco to the prefix (turned on
       # by --with-draco without argument)
       if test "${with_draco}" = yes ; then
	   
	   # set DRACO_INC and DRACO_LIB
	   DRACO_INC="${prefix}/include"
	   DRACO_LIB="${prefix}/lib"

       # otherwise it is the directory where the installed lib and 
       # include are; check them
       else
	   
	   # set DRACO_INC and DRACO_LIB
	   DRACO_INC="${with_draco}/include"
	   DRACO_LIB="${with_draco}/lib"

	   # make sure they exist
	   if test ! -d "${DRACO_INC}" ; then
	       AS_MSG_ERROR("${DRACO_INC} does not exist")
	   elif test ! -d "${DRACO_LIB}" ; then
	       AC_MSG_ERROR("${DRACO_LIB} does not exist")
	   fi

       fi

   fi

   # check to see if --with_draco_inc or --with_draco_lib have 
   # been set
   else
     
       # if --with_draco_inc is defined then set and check DRACO_INC
       if test -n "${with_draco_inc}" ; then
	   DRACO_INC="${with_draco_inc}"
	   if test ! -d "${DRACO_INC}" ; then
	       AC_MSG_ERROR("${DRACO_INC} does not exist")
	   fi
       fi

       # if --with_draco_lib is defined then set and check DRACO_LIB
       if test -n "${with_draco_lib}" ; then
	   DRACO_LIB="${with_draco_lib}"
	   if test ! -d "${DRACO_LIB}" ; then
	       AC_MSG_ERROR("${DRACO_LIB} does not exist")
	   fi
       fi

   fi

   # check settings of DRACO_INC and DRACO_LIB
   if test -n "${DRACO_INC}" && test -n "${DRACO_LIB}" ; then
       AC_MSG_RESULT{"${DRACO_INC} and ${DRACO_LIB} set")
   elif test -n "${DRACO_INC}" && test -z "${DRACO_LIB}" ; then
       AC_MSG_RESULT("${DRACO_INC} set, draco libraries will be searched in standard locations")
   elif test -n "${DRACO_LIB}" && test -z "${DRACO_INC}" ; then
       AC_MSG_RESULT("${DRACO_LIB} set, draco includes will be searched in standard locations")
   else
       AC_MSG_RESULT("use default search paths for draco libraries and includes")

   # add draco include directory to CPPFLAGS
   if test -n "${DRACO_INC}" ; then
       CPPFLAGS="${CPPFLAGS} -I ${DRACO_INC}"
   fi
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_draco_vendor.m4
dnl-------------------------------------------------------------------------dnl
