dnl-------------------------------------------------------------------------dnl
dnl ac_local.m4
dnl
dnl Macros used internally within the Draco build system.
dnl
dnl Thomas M. Evans
dnl 1999/02/04 01:56:22
dnl-------------------------------------------------------------------------dnl
##---------------------------------------------------------------------------##
## $Id$
##---------------------------------------------------------------------------##

dnl KT (2010-04-27): Replace these bits of logic with something that
dnl looks more like this?

dnl AC_ARG_WITH(xmlparser,
dnl     AS_HELP_STRING([--with-xmlparser], [select used xml parser; you can chose from expat,
dnl 	libxml2 or auto. @<:@default=libxml2@:>@]),
dnl     [case "$withval" in
dnl         expat)    CONFIG_XMLPARSER=expat ;;
dnl         libxml2)  CONFIG_XMLPARSER=libxml2 ;;
dnl         n | no)   CONFIG_XMLPARSER=no ;;
dnl         auto | *) CONFIG_XMLPARSER=auto ;;
dnl     esac],
dnl     [CONFIG_XMLPARSER=auto]
dnl )

dnl-------------------------------------------------------------------------dnl
dnl AC_SETUP_VENDOR
dnl
dnl Defines three command line options:
dnl    --with-xxx=[no|default_value]
dnl    --with-xxx-inc=[DIR]
dnl    --with-xxx-lib=[DIR]
dnl 
dnl If no command line option is provided, value defaults to 'no'
dnl (equivalent to --without-xxx).  If option is provided without
dnl arguments, then the value 'yes' is assigned temporarily.  The
dnl following values will be used in this case:
dnl    with_xxx     <-- $3
dnl    with_xxx_inc <-- $XXX_INC_DIR
dnl    with_xxx_lib <-- $XXX_LIB_DIR
dnl
dnl Usage: AC_SETUP_VENDOR( vendor_name, no, default_value,
dnl                        @S|@XXX_INC_DIR, @S|@XXX_LIB_DIR )
dnl 
dnl     $1 is the vendor name
dnl     $2 default ${with_vendor} to on (yes) or off (no) if user
dnl        doesn't specify a value.
dnl     $3 is the library name (e.g.: grace_np)
dnl     $4 is a string for the help message.
dnl     $5 is a string for the help message.
dnl     $6 is a string for the help message (known valid option values)
dnl
dnl Example: AC_SETUP_VENDOR( grace, no, grace_np, 
dnl                          @S|@{GRACE_INC_DIR}, @S|@{GRACE_LIB_DIR} )
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_SETUP_VENDOR], [dnl

   dnl echo "in setup_vendor: $1"

   dnl Vendor name with underscores in place of hyphens
   dnl creates '--with-xxx' option
   define([AC_VENDOR], [translit($1, [-], [_])])dnl
   dnl This name used to generate '--with-xxx-inc'
   define([AC_VENDOR_INC], AC_VENDOR[_inc])dnl
   dnl This name used to generate '--with-xxx-lib'
   define([AC_VENDOR_LIB], AC_VENDOR[_lib])dnl
   dnl Define the variable name ${with_xxx}
   define([AC_WITH_VENDOR], [with_]AC_VENDOR)dnl
   dnl Define the variable name ${with_xxx_inc}
   define([AC_WITH_VENDOR_INC], [with_]AC_VENDOR_INC)dnl
   dnl Define the variable name ${with_xxx_lib}
   define([AC_WITH_VENDOR_LIB], [with_]AC_VENDOR_LIB)dnl
   dnl Define the environment variable name for the include path, 
   dnl ${XXX_INC_DIR}
   pushdef([VENDOR_INC_DIR], patsubst(translit([$1_inc_dir], [a-z], [A-Z]), -, _))dnl
   dnl Define the environment variable name for the library path, 
   dnl ${XXX_LIB_DIR}
   pushdef([VENDOR_LIB_DIR], patsubst(translit([$1_lib_dir], [a-z], [A-Z]), -, _))dnl
   define([AC_CMDLINE],dnl
      [echo "$]AC_WITH_VENDOR[" | sed 's%//*%/%g' | sed 's%/$%%'])dnl
   define([AC_CMDLINE_INC],dnl
      [echo "$]AC_WITH_VENDOR_INC[" | sed 's%//*%/%g' | sed 's%/$%%'])dnl
   define([AC_CMDLINE_LIB],dnl
      [echo "$]AC_WITH_VENDOR_LIB[" | sed 's%//*%/%g' | sed 's%/$%%'])dnl

   dnl ------------------------------------------------------------ dnl

   # Keep track of user provided commands. 
   vendor_requested_by_user=false

   dnl define --with-<vendor> option. Produces the help message:
   dnl --with-<vendor>[=[no|libname]]
   dnl             use vendor <vendor> (libname) library (disabled by default)
   AC_ARG_WITH([$1],
   [AS_HELP_STRING( [--with-$1@<:@=@<:@$2|$6@:>@@:>@],
       [use vendor $1 ($3) library (disabled by default)])],

   dnl Execute this block if --with-xxx is provided. 
   vendor_requested_by_user=true
   dnl   echo "   setting vendor_requested_by_user to true"
   if test $AC_WITH_VENDOR != "no" ; then
      if test $AC_WITH_VENDOR = "yes" ; then
         # following eval needed to remove possible '\' from $3
         eval AC_WITH_VENDOR=$3
      fi
      # this command removes double slashes and any trailing slash
      AC_WITH_VENDOR=`eval AC_CMDLINE`
      if test "$AC_WITH_VENDOR:-null}" = "null" ; then
         { echo "configure: error: --with-$1 directory is unset" 1>&2; \
           exit 1; }
      fi
      # this sets up the shell variable, with the name of the CPPtoken,
      # and that we later will do an AC_SUBST on.
      dnl kt $2="${AC_WITH_VENDOR}/"
      dnl AC_WITH_VENDOR="${AC_WITH_VENDOR}/"
      dnl this defines the CPP macro with the directory and single slash appended.
      dnl kt AC_DEFINE_UNQUOTED($2, ${AC_WITH_VENDOR}/)dnl
      dnl print a message to the users (that can be turned off with --silent)
      dnl  echo "AC_WITH_VENDOR has been set to $AC_WITH_VENDOR" 1>&6
   fi,

   dnl Execute this block if --with-xxx is not provided. 
      if test "$2" = "yes"; then
        AC_WITH_VENDOR=$3
      else
        AC_WITH_VENDOR=no
      fi
   )

   dnl ------------------------------------------------------------ dnl

   dnl define --with-<vendor>-inc:

   dnl If '--with-<vendor>' is not provided, but either or both
   dnl '--with-<vendor>-lib' and '--with-<vendor>-inc' are provided,
   dnl then enable vendor_<vendor>.

   if ! test "$4x" = "x"; then

   AC_ARG_WITH([$1-inc],
   [AS_HELP_STRING( [--with-$1-inc@<:@=@<:@no|DIR@:>@@:>@],
       [specify the location of the $1 header files (disabled by ]
       [default, if yes, $4 will be used.)])], 

   dnl Execute this block if --with-xxx-inc is provided.
   if test $AC_WITH_VENDOR_INC != "no" ; then
      # Turn on current vendor unless user assigned "no" and assign
      # library name  
      if test AC_WITH_VENDOR = "no"; then
         { echo "configure: error: option --without-$1 conflicts with option --with-$1-inc" 1>&2; \
           exit 1; }
      else
         eval AC_WITH_VENDOR=$3
      fi
      # If user turns on option, but does not provide a path, use
      # value provided by ${XXX_INC_DIR}
      if test $AC_WITH_VENDOR_INC = "yes" ; then
         # following eval needed to remove possible '\' from the
         # environment variriable ${XXX_INC_DIR}
         eval AC_WITH_VENDOR_INC=$VENDOR_INC_DIR
      fi
      # this command removes double slashes and any trailing slash
      AC_WITH_VENDOR_INC=`eval AC_CMDLINE_INC`
      # Error if no value available for with_xxx_inc
      if test "${AC_WITH_VENDOR_INC:-null}" = "null" ; then
         { echo "configure: error: --with-$1-inc directory is unset" 1>&2; \
           exit 1; }
      fi
      # Error if value provided by with_xxx_inc is not valid.
      if test ! -d ${AC_WITH_VENDOR_INC} ; then
         { echo "configure: error: $AC_WITH_VENDOR_INC: invalid directory" 1>&2; \
           exit 1; }
      fi
      dnl print a message to the users (that can be turned off with --silent)
      dnl echo "AC_WITH_VENDOR_INC has been set to $AC_WITH_VENDOR_INC" 1>&6
   fi,

   dnl Execute this block if --with-xxx-inc is not provided.
     dnl echo "Using default (no) for AC_WITH_VENDOR_INC"
     AC_WITH_VENDOR_INC=no
     dnl If the value for with_<vendor>_inc is anything but 'no' (the
     dnl deafult) then turn on <vendor> and try to figure out the
     dnl include path.
     if test ${vendor_requested_by_user} = false || ! test $AC_WITH_VENDOR = "no"; then
        eval AC_WITH_VENDOR_INC=${VENDOR_INC_DIR}
       dnl echo "   Setting AC_WITH_VENDOR_INC to VENDOR_INC_DIR = ${VENDOR_INC_DIR}"
        # Turn off vendor if no value provided or if an invalid directory
        # is provided..
        if test "${AC_WITH_VENDOR_INC:-null}" = "null" || 
           test ! -d ${AC_WITH_VENDOR_INC}; then
           AC_WITH_VENDOR_INC=no
           AC_WITH_VENDOR=no
           if ${vendor_requested_by_user} = true; then
             { echo "configure: error: --with-$1 requested but vendor not found (inc)." 1>&2; \
              exit 1; }
           fi
        fi
     fi

   )
   fi

   dnl ------------------------------------------------------------ dnl


   dnl define --with-<vendor>-lib:

   dnl If '--with-<vendor>' is not provided, but either or both
   dnl '--with-<vendor>-lib' and '--with-<vendor>-lib' are provided,
   dnl then enable vendor_<vendor>.
   AC_ARG_WITH([$1-lib],
   [AS_HELP_STRING( [--with-$1-lib@<:@=@<:@no|DIR@:>@@:>@],
       [specify the location of the $1 header files (disabled by ]
       [default, if yes, $5 will be used.)])],

   dnl Execute this block if --with-xxx-lib is provided.
   if test $AC_WITH_VENDOR_LIB != "no" ; then
      # Turn on current vendor unless user assigned "no" and assign
      # library name  
      if test AC_WITH_VENDOR = "no"; then
         { echo "configure: error: option --without-$1 conflicts with option --with-$1-lib" 1>&2; \
           exit 1; }
      else
         eval AC_WITH_VENDOR=$3
      fi
      # If user turns on option, but does not provide a path, use
      # value provided by ${XXX_LIb_DIR}
      if test $AC_WITH_VENDOR_LIB = "yes" ; then
         # following eval needed to remove possible '\' from the
         # environment variriable ${XXX_LIB_DIR}
         eval AC_WITH_VENDOR_LIB=$VENDOR_LIB_DIR
      fi
      # this command removes double slashes and any trailing slash
      AC_WITH_VENDOR_LIB=`eval AC_CMDLINE_LIB`
      # Error if no value available for with_xxx_lib
      if test "${AC_WITH_VENDOR_LIB:-null}" = "null" ; then
         { echo "configure: error: --with-$1-lib directory is unset" 1>&2; \
           exit 1; }
      fi
      # Error if value provided by with_xxx_lib is not valid.
      if test ! -d ${AC_WITH_VENDOR_LIB} ; then
         { echo "configure: error: ${AC_WITH_VENDOR_LIB}: invalid directory" 1>&2; \
           exit 1; }
      fi
   fi,

   dnl Execute this block if --with-xxx-lib is not provided.
   dnl When no command line value is provided
     dnl echo "   Using default (no) for AC_WITH_VENDOR_LIB"
     AC_WITH_VENDOR_LIB=no
     dnl If the value for with_<vendor>_lib is anything but 'no' (the
     dnl deafult) then turn on <vendor> and try to figure out the
     dnl include path.
     if test ${vendor_requested_by_user} = false || ! test $AC_WITH_VENDOR = "no"; then
        eval AC_WITH_VENDOR_LIB=${VENDOR_LIB_DIR}
        dnl echo "   Setting AC_WITH_VENDOR_LIB to VENDOR_LIB_DIR = ${VENDOR_LIB_DIR}"
        # Turn off vendor if no value provided or if an invalid directory
        # is provided..
        if test "${AC_WITH_VENDOR_LIB:-null}" = "null" || 
           test ! -d ${AC_WITH_VENDOR_LIB}; then
           dnl echo "   Bad library path. Turning off vendor."
           AC_WITH_VENDOR_LIB=no
           AC_WITH_VENDOR=no
           if ${vendor_requested_by_user} = true; then
             { echo "configure: error: --with-$1 requested but vendor not found (lib)." 1>&2; \
              exit 1; }
           fi
        fi
     fi
   )

    dnl echo "   setup_vendor: VENDOR      = AC_VENDOR"
    dnl echo "   setup_vendor: WITH_VENDOR = AC_WITH_VENDOR = ${AC_WITH_VENDOR}"
    dnl echo "   setup_vendor: WITH_VENDOR_INC = AC_WITH_VENDOR_INC = ${AC_WITH_VENDOR_INC}"
    dnl echo "   setup_vendor: WITH_VENDOR_LIB = AC_WITH_VENDOR_LIB = ${AC_WITH_VENDOR_LIB}"

])

dnl-------------------------------------------------------------------------dnl
dnl AC_FINALIZE_VENDOR
dnl
dnl Once the vendor has been found and verified, update the build
dnl system variables: VENDOR_INC, VENDOR_LIBS, VENDOR_TEST_LIBS,
dnl VENDOR_LIB_DIRS and VENDOR_INC_DIRS.
dnl
dnl Usage: AC_FINALIZE_VENDOR( vendor_name, vendor_extra_libs )
dnl
dnl     $1 is the vendor name (e.g.: grace or gsl)
dnl
dnl Example: AC_FINALIZE_VENDOR( grace )
dnl          AC_FINALIZE_VENDOR( [gsl], [${gsl_extra_libs}] )
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_FINALIZE_VENDOR], [dnl

   define([AC_VENDOR], [translit(m4_strip($1),[-],[_])])dnl     e.g.: grace
   define([AC_WITH_VENDOR], [with_]AC_VENDOR)dnl      e.g.: with_grace
   define([AC_VENDOR_VENDOR], [vendor_]AC_VENDOR)dnl  e.g.: vendor_grace
   define([AC_VENDOR_INC], [with_]AC_VENDOR[_inc])dnl e.g.: with_grace_inc
   define([AC_VENDOR_LIB], [with_]AC_VENDOR[_lib])dnl e.g.: with_grace_lib
   pushdef([VENDOR_LIB_DIR], translit(m4_strip($1)[_lib_dir], [a-z], [A-Z]))dnl

   dnl echo "in finalize_vendor: (m4_strip($1))"
   dnl echo "   VENDOR_VENDOR = AC_VENDOR_VENDOR = ${AC_VENDOR_VENDOR}"

   vendor_extra_libs=$2

   # set up the libraries and include path
   if test -n "${AC_VENDOR_VENDOR}" ; then
       dnl echo "VENDOR_INC  = AC_VENDOR_INC  = ${AC_VENDOR_INC}"
       dnl echo "VENDOR_LIB  = AC_VENDOR_LIB  = ${AC_VENDOR_LIB}"
       dnl echo "WITH_VENDOR = AC_WITH_VENDOR = ${AC_WITH_VENDOR}"
       
       # include path
       if test -d "${AC_VENDOR_INC}"; then
           # add to include path
           dnl echo "   AC_VENDOR_INC valid, adding..."
           VENDOR_INC="${VENDOR_INC} -I${AC_VENDOR_INC}"
       fi

       # library path
       if test -d "${AC_VENDOR_LIB}"; then
          dnl echo "   AC_VENDOR_LIB valid, adding..."
          if test "${AC_VENDOR_VENDOR}" = pkg; then
             dnl echo "   pkg -L -l"
             # if with_vendor=='vendor', don't add anything (should be automatic)
             if ! test "${AC_WITH_VENDOR}" = "vendor"; then
                VENDOR_LIBS="${VENDOR_LIBS} -L${AC_VENDOR_LIB} -l${AC_WITH_VENDOR}"
             fi
             if ! test "${vendor_extra_libs}x" = "x"; then
                VENDOR_LIBS="${VENDOR_LIBS} ${vendor_extra_libs}"
             fi
          elif test "${AC_VENDOR_VENDOR}" = test; then
             dnl echo "   test -L -l"
             # if with_vendor=='vendor', don't add anything (should be automatic)
             if ! test "${AC_WITH_VENDOR}" = "vendor"; then
                VENDOR_TEST_LIBS="${VENDOR_TEST_LIBS} -L${AC_VENDOR_LIB} -l${AC_WITH_VENDOR}"
             fi
             if ! test "${vendor_extra_libs}x" = "x"; then
                VENDOR_TEST_LIBS="${VENDOR_TEST_LIBS} ${vendor_extra_libs}"
             fi
          fi
       else
          { echo "configure: error: --with-m4_strip($1)=${AC_WITH_VENDOR} requested but vendor " 1>&2; \
            echo "   library directory not provided. You must provide the directory location by " 1>&2; \
	    echo "   setting VENDOR_LIB_DIR in your environment or by providing the path on " 1>&2; \
 	    echo "   the configure command line with the option --with-m4_strip($1)-lib=<path>." 1>&2; \
            exit 1; }
dnl                echo "failed link test.  Try setting --with-$1-lib=PATH or setting" 1>&2; \
dnl                echo "@S|@AC_VENDOR_LIB." 1>&2; \
dnl                exit 1; }
dnl          echo "AC_VENDOR_LIB not valid, only add -l...($7)"
          dnl check to see if we can find the requested library
          dnl vendor_lib_found=false
	  dnl mylib="${AC_WITH_VENDOR}"
dnl 	  myfunction="foo"
dnl           if test -n $7; then
dnl             myfunction=$7
dnl           fi
dnl 	  AC_CHECK_LIB( [${mylib}], [${myfunction}], [vendor_lib_found=true] )
dnl           if test "${vendor_lib_found}" = true; then
dnl              if test "${AC_VENDOR_VENDOR}" = pkg; then
dnl                 echo "pkg -l"
dnl                 VENDOR_LIBS="${VENDOR_LIBS} -l${AC_WITH_VENDOR}"
dnl              elif test "${AC_VENDOR_VENDOR}" = test; then
dnl                 echo "test -l"
dnl                 VENDOR_TEST_LIBS="${VENDOR_TEST_LIBS} -l${AC_WITH_VENDOR}"
dnl              fi
dnl           else
dnl              { echo "configure: error: --with-$1 requested but vendor library " 1>&2; \
dnl                echo "failed link test.  Try setting --with-$1-lib=PATH or setting" 1>&2; \
dnl                echo "@S|@AC_VENDOR_LIB." 1>&2; \
dnl                exit 1; }

       fi

       # add VENDOR directory to VENDOR_LIB_DIRS
       VENDOR_LIB_DIRS="${VENDOR_LIB_DIRS} ${AC_VENDOR_LIB}"
       VENDOR_INC_DIRS="${VENDOR_INC_DIRS} ${AC_VENDOR_INC}"

       dnl echo "   VENDOR_INC       = ${VENDOR_INC}"
       dnl echo "   VENDOR_LIBS      = ${VENDOR_LIBS}"
       dnl echo "   VENDOR_TEST_LIBS = ${VENDOR_TEST_LIBS}"

   fi


])
dnl-------------------------------------------------------------------------dnl
dnl AC_NO_VENDOR_WARN
dnl
dnl Print a warning message if a package cannot be configured due to a
dnl missing vendor package.
dnl
dnl AC_NO_VENDOR_WARN( package_name, vendor_name )
dnl
dnl used by: src/configure.ac
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_NO_VENDOR_WARN], [dnl

   define([AC_PACKAGE_NAME],[$1])
   define([AC_VENDOR_NAME],[$2])
dnl   pushdef([AC_INC_DIR], translit([$2[_inc_dir]], [a-z], [A-Z]))dnl
dnl   pushdef([AC_LIB_DIR], translit([$2[_lib_dir]], [a-z], [A-Z]))dnl

   AC_MSG_WARN([Package AC_PACKAGE_NAME will not be built because \
vendor AC_VENDOR_NAME was excluded or not found.  Try configure options \
--with-AC_VENDOR_NAME=<lib>, --with-AC_VENDOR_NAME-lib=DIR and --with-AC_VENDOR_NAME-inc=DIR. \
You may also want to set the environment variables @S|@<VENDOR>_INC_DIR \
and @S|@<VENDOR>_LIB_DIR.]) 
])

dnl-------------------------------------------------------------------------dnl
dnl AC_WITH_DIR
dnl
dnl Define --with-xxx[=DIR] with defaults to an environment variable.
dnl       Usage: AC_WITH_DIR(flag, CPPtoken, DefaultValue, HelpStr)
dnl                for environment variables enter \${ENVIRONVAR} for
dnl                DefaultValue
dnl usage: in aclocal.m4
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_WITH_DIR], [dnl

 dnl
 dnl The following M4 macros will be expanded into the body of AC_ARG_WITH
 dnl
 dnl AC_PACKAGE is the flag with all dashes turned to underscores
 dnl AC_WITH_PACKAGE will be substituted to the autoconf shell variable
 dnl    with_xxx
 dnl AC_CMDLINE is the shell command to strip double and trailing slashes
 dnl    from directory names.

 define([AC_PACKAGE], [translit($1, [-], [_])])dnl
 define([AC_WITH_PACKAGE], [with_]AC_PACKAGE)dnl
 define([AC_CMDLINE],dnl
[echo "$]AC_WITH_PACKAGE[" | sed 's%//*%/%g' | sed 's%/$%%'])dnl

dnl   [  --with-$1[=DIR]    $4 ($3 by default)],
 AC_ARG_WITH([$1],
   [AS_HELP_STRING([--with-$1@<:@=DIR@:>@],[$4 ($3 by default)])],
   if test $AC_WITH_PACKAGE != "no" ; then
      if test $AC_WITH_PACKAGE = "yes" ; then
         # following eval needed to remove possible '\' from $3
         eval AC_WITH_PACKAGE=$3
      fi

      # this command removes double slashes and any trailing slash

      AC_WITH_PACKAGE=`eval AC_CMDLINE`
      if test "$AC_WITH_PACKAGE:-null}" = "null" ; then
         { echo "configure: error: --with-$1 directory is unset" 1>&2; \
           exit 1; }
      fi
      if test ! -d $AC_WITH_PACKAGE ; then
         { echo "configure: error: $AC_WITH_PACKAGE: invalid directory" 1>&2; \
           exit 1; }
      fi

      # this sets up the shell variable, with the name of the CPPtoken,
      # and that we later will do an AC_SUBST on.
      $2="${AC_WITH_PACKAGE}/"

      # this defines the CPP macro with the directory and single slash appended.
      AC_DEFINE_UNQUOTED($2, ${AC_WITH_PACKAGE}/)dnl

      # print a message to the users (that can be turned off with --silent)

      dnl echo "$2 has been set to $$2" 1>&6

   fi)

   AC_SUBST($2)dnl

])
	
dnl-------------------------------------------------------------------------dnl
dnl AC_VENDORLIB_SETUP(1,2)
dnl
dnl set up for VENDOR_LIBS or VENDOR_TEST_LIBS
dnl usage: in aclocal.m4
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_VENDORLIB_SETUP], [dnl

   # $1 is the vendor_<> tag (equals pkg or test)
   # $2 are the directories added 

   if test "${$1}" = pkg ; then
       VENDOR_LIBS="${VENDOR_LIBS} $2"
   elif test "${$1}" = test ; then
       VENDOR_TEST_LIBS="${VENDOR_TEST_LIBS} $2"
   fi
])

dnl-------------------------------------------------------------------------dnl
dnl AC_FIND_TOP_SRC(1,2)
dnl 
dnl Find the top source directory of the package by searching upward
dnl from the argument directory. The top source directory is defined
dnl as the one with a 'config' sub-directory.
dnl
dnl Note: This function will eventually quit if the searched for
dnl directory is not above the argument. It does so when $temp_dir
dnl ceases to be a valid directory, which only seems to happen after a
dnl LOT of ..'s are added to it.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_FIND_TOP_SRC], [dnl
   
   # $1 is the component's source directory
   # $2 is the variable to store the package's main source directory in.

   temp_dir=$1
   AC_MSG_CHECKING([package top source directory])
   while test -d $temp_dir -a ! -d $temp_dir/config ; do   
       temp_dir="${temp_dir}/.."
   done
   if test -d $temp_dir; then
       $2=`cd $temp_dir; pwd;`
       AC_MSG_RESULT([$$2])
   else
       AC_MSG_ERROR('Could not find package top source directory')
   fi
])

dnl-------------------------------------------------------------------------dnl
dnl DO VARIABLE SUBSTITUTIONS ON AC_OUTPUT
dnl
dnl These are all the variable substitutions used within the draco
dnl build system
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DBS_VAR_SUBSTITUTIONS], [dnl

   # these variables are declared "precious", meaning that they are
   # automatically substituted, put in the configure --help, and
   # cached 
   AC_ARG_VAR(CC)dnl
   AC_ARG_VAR(CFLAGS)dnl

   AC_ARG_VAR(CXX)dnl
   AC_ARG_VAR(CXXFLAGS)dnl

   AC_ARG_VAR(LD)dnl
   AC_ARG_VAR(LDFLAGS)dnl

   AC_ARG_VAR(AR)dnl
   AC_ARG_VAR(ARFLAGS)dnl

   AC_ARG_VAR(CPPFLAGS)dnl

   # dependency rules
   AC_SUBST(DEPENDENCY_RULES)

   # other compiler substitutions
   AC_SUBST(STRICTFLAG)dnl
   AC_SUBST(PARALLEL_FLAG)dnl
   AC_SUBST(RPATH)dnl
   AC_SUBST(LIB_PREFIX)dnl

   # install program
   AC_SUBST(INSTALL)dnl
   AC_SUBST(INSTALL_DATA)dnl

   # files to install
   : ${installfiles:='${install_executable} ${install_lib} ${install_headers}'}
   AC_SUBST(installfiles)dnl
   AC_SUBST(install_executable)dnl
   AC_SUBST(install_lib)dnl
   AC_SUBST(install_headers)dnl
   AC_SUBST(installdirs)dnl

   # package libraries
   AC_SUBST(alltarget)dnl
   AC_SUBST(libsuffix)dnl
   AC_SUBST(dirstoclean)dnl
   AC_SUBST(package)dnl
   AC_SUBST(DRACO_DEPENDS)dnl
   AC_SUBST(DRACO_LIBS)dnl
   AC_SUBST(VENDOR_DEPENDS)dnl
   AC_SUBST(VENDOR_INC)dnl
   AC_SUBST(VENDOR_LIBS)dnl
   AC_SUBST(ARLIBS)dnl

   # package testing libraries
   AC_SUBST(PKG_DEPENDS)dnl
   AC_SUBST(PKG_LIBS)dnl
   AC_SUBST(DRACO_TEST_DEPENDS)dnl
   AC_SUBST(DRACO_TEST_LIBS)dnl
   AC_SUBST(VENDOR_TEST_DEPENDS)dnl
   AC_SUBST(VENDOR_TEST_LIBS)dnl
   AC_SUBST(ARTESTLIBS)dnl
   AC_SUBST(test_flags)dnl
   AC_SUBST(test_scalar)dnl
   AC_SUBST(test_nprocs)dnl
   AC_SUBST(test_output_files)dnl
   AC_SUBST(scalar_tests)dnl
   AC_SUBST(parallel_tests)dnl
   AC_SUBST(app_tests)dnl
   AC_SUBST(app_test_nprocs)dnl

   # libraries
   AC_ARG_VAR(LIBS)dnl

   # configure options
   AC_SUBST(configure_command)dnl

   # directories in source tree
   AC_SUBST(package_top_srcdir)
   
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_local.m4
dnl-------------------------------------------------------------------------dnl
