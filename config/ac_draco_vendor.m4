dnl-------------------------------------------------------------------------dnl
dnl ac_draco_vendor.m4
dnl macros for setting up draco as a vendor
dnl
dnl $Id$
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl AC_DRACO_VENDOR_SETUP
dnl
dnl Setup macro that allows draco to be defined as a vendor.
dnl This should be called in a draco-client configure.in as follows:
dnl
dnl AC_DRACO_VENDOR_SETUP
dnl
dnl AC_NEEDS_LIBS(pkg_components)
dnl AC_NEEDS_LIBS_DRACO(cdi ds++)
dnl
dnl That is, it must be called before AC_NEEDS_LIBS_DRACO.
dnl
dnl Three options are defined that tell configure where draco includes
dnl and libraries are located:
dnl
dnl  --with-draco
dnl  --with-draco-lib
dnl  --with-draco-inc
dnl
dnl Setting --with-draco=<dir> tells configure to look in <dir>/lib,
dnl <dir>/inc, <dir>/libexec for draco stuff.  The include and/or lib 
dnl directories can be overridden by setting --with-draco-inc or
dnl --with-draco-lib.  If --with-draco is set, without arguments, then
dnl --prefix is used as the draco directory.  The default location (if
dnl --with-draco is not called at all) is /usr/local.  Of course, this
dnl can be overwritten for the lib and include directories by using 
dnl --with-draco-lib and/or --with-draco-inc (libexec is found in the
dnl directory set by --with-draco-lib (in libexec instead of lib.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DRACO_VENDOR_SETUP, [dnl

   dnl define --with-draco 
   AC_ARG_WITH(draco,
      [  --with-draco            give head location of draco includes and libs (default is /usr/local; no argmument sets to prefix)])

   AC_ARG_WITH(draco-inc,
      [  --with-draco-inc        give location of draco includes])

   AC_ARG_WITH(draco-lib,
      [  --with-draco-lib        give location of draco libs])

   # draco include and lib dirs
   DRACO_INC=''
   DRACO_LIB=''
   draco_in_prefix=''

   AC_MSG_CHECKING("draco directories") 

   # if we use the with-draco option to specify the location of draco,
   # then check it
   if test -n "${with_draco}" ; then 

       # check to see if we should set draco to the prefix (turned on
       # by --with-draco without argument)
       if test "${with_draco}" = yes ; then
	   
	   # set DRACO_INC and DRACO_LIB
	   DRACO_INC="${prefix}/include"
	   DRACO_LIB="${prefix}/lib"
	   libexecdir="${prefix}/libexec"

	   draco_in_prefix='true'

       # otherwise it is the directory where the installed lib and 
       # include are; check them
       else
	   
	   # set DRACO_INC and DRACO_LIB
	   DRACO_LIB="${with_draco}/lib"   
	   DRACO_INC="${with_draco}/include"
	   libexecdir="${with_draco}/libexec"

       fi

   # set draco default location to /usr/local
   else

       DRACO_INC="/usr/local/include"
       DRACO_LIB="/usr/local/lib"
       libexecdir="/usr/local/libexec"

   fi
     
   # if --with_draco_inc is defined then set and check DRACO_INC
   if test -n "${with_draco_inc}" ; then
       DRACO_INC="${with_draco_inc}"
   fi

   # if --with_draco_lib is defined then set and check DRACO_LIB
   if test -n "${with_draco_lib}" ; then
       DRACO_LIB="${with_draco_lib}"
       libexecdir="${with_draco_lib}/../libexec"
   fi

   # make sure they exist   
   if test ! -d "${DRACO_INC}" && test -z "${draco_in_prefix}" ; then
       AC_MSG_ERROR("${DRACO_INC} does not exist")
   fi

   if test ! -d "${DRACO_LIB}" && test -z "${draco_in_prefix}" ; then
       AC_MSG_ERROR("${DRACO_LIB} does not exist")
   fi


   AC_MSG_RESULT("${DRACO_INC} and ${DRACO_LIB} and ${libexecdir} set")

   # add draco include directory to CPPFLAGS
   if test -z "${draco_in_prefix}" ; then
       CPPFLAGS="${CPPFLAGS} -I${DRACO_INC}"
   fi

   # add draco to VENDIR_DIRS
   VENDOR_DIRS="${DRACO_LIB} ${VENDOR_DIRS}"
])

dnl-------------------------------------------------------------------------dnl
dnl AC_NEEDS_LIBS_DRACO
dnl
dnl add DRACO libraries necessary for a package
dnl this MACRO must be called after AC_NEEDS_LIBS
dnl usage: configure.in
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_NEEDS_LIBS_DRACO, [dnl
   if test ${has_draco_libdir:=no} != "yes" ; then
       
       has_draco_libdir="yes"

       if test ${has_libdir} != "yes" ||
	  test ${draco_in_prefix} != "true" ; then
	       DRACO_LIBS="${DRACO_LIBS} -L${DRACO_LIB}"
       fi

   fi

   for lib in $1
   do
       # temporary string to keep line from getting too long
       draco_depends="${DRACO_LIB}/lib${lib}\${libsuffix}"
       DRACO_DEPENDS="${DRACO_DEPENDS} ${draco_depends}"
       DRACO_LIBS="${DRACO_LIBS} -l${lib}"
   done
])

dnl-------------------------------------------------------------------------dnl
dnl AC_NEEDS_LIBS_TEST_DRACO
dnl
dnl add DRACO libraries necessary for testing a package
dnl this MACRO must be called after AC_NEEDS_LIBS_TEST
dnl usage: configure.in
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_NEEDS_LIBS_TEST_DRACO, [dnl
   DRACO_TEST_LIBS="${DRACO_TEST_LIBS} -L${DRACO_LIB}"

   for lib in $1
   do
       # temporary string to keep line from getting too long
       draco_test_depends="${DRACO_LIB}/lib${lib}\${libsuffix}"
       DRACO_TEST_DEPENDS="${DRACO_TEST_DEPENDS} ${draco_test_depends}"
       DRACO_TEST_LIBS="${DRACO_TEST_LIBS} -l${lib}"
   done
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_draco_vendor.m4
dnl-------------------------------------------------------------------------dnl
