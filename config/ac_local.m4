dnl-------------------------------------------------------------------------dnl
dnl ac_local.m4
dnl service macros used in ac_vendors.m4, ac_dracoarg.m4, and ac_dracoenv.m4
dnl
dnl $Id$
dnl-------------------------------------------------------------------------dnl

dnl-------------------------------------------------------------------------dnl
dnl AC_WITH_DIR
dnl
dnl Define --with-xxx[=DIR] with defaults to an environment variable.
dnl       Usage: AC_WITH_DIR(flag, CPPtoken, DefaultValue, HelpStr)
dnl                for environment variables enter \${ENVIRONVAR} for
dnl                DefaultValue
dnl usage: in aclocal.m4
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_WITH_DIR, [dnl

 dnl
 dnl  The following M4 macros will be expanded into the body of AC_ARG_WITH
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

 AC_ARG_WITH($1,
   [  --with-$1[=DIR]    $4 ($3 by default)],
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

      echo "$2 has been set to $$2" 1>&6

   fi)

   AC_SUBST($2)dnl

])
	
dnl-------------------------------------------------------------------------dnl
dnl AC_VENDORLIB_SETUP(1,2)
dnl
dnl set up for VENDOR_LIBS or VENDOR_TEST_LIBS
dnl usage: in aclocal.m4
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_VENDORLIB_SETUP, [dnl

   # $1 is the vendor_<> tag (equals pkg or test)
   # $2 are the directories added 

   if test "${$1}" = pkg ; then
       VENDOR_LIBS="${VENDOR_LIBS} $2"
   elif test "${$1}" = test ; then
       VENDOR_TEST_LIBS="${VENDOR_TEST_LIBS} $2"
   fi
])

dnl-------------------------------------------------------------------------dnl
dnl AC_DETERMINE_INT
dnl
dnl DETERMINE C++ DATA TYPE FOR A GIVEN INTEGER SIZE
dnl eg. AC_DETERMINE_INT(4) sets the variable INTEGER_SIZE_TYPE to int
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DETERMINE_INT, [dnl

   AC_MSG_CHECKING("C++ data type of integer of size $1 bytes")

   # set the language to C++
   AC_LANG_CPLUSPLUS

   check_ints='false'
   
   # check to see if regular int does the job
   AC_TRY_RUN([
       int main()
       {
	   int p = 1;
	   if (sizeof(int) == $1)
	       p = 0;
	   return p;
       }], [INTEGER_SIZE_TYPE='int'], 
	   [check_ints='true'])
      
   # check long int type if ints did not pass
   if test "${check_ints}" = true ; then

       AC_TRY_RUN([
	   int main()
	   {   
	       int p = 1;
	       if (sizeof(long) == $1)
		   p = 0;
	       return p;
	   }], [check_ints='false';
		INTEGER_SIZE_TYPE='long'], [])

   fi  
     
   # check long long type if long did not pass
   if test "${check_ints}" = true ; then

       AC_TRY_RUN([
	   int main()
	   {   
	       int p = 1;
	       if (sizeof(long long) == $1)
		   p = 0;
	       return p;
	   }], [check_ints='false';
		INTEGER_SIZE_TYPE='long long'], [])

   fi
       
   # error message because we haven't found a valid type
   if test "${check_ints}" = true ; then
       AC_MSG_RESULT("no match found")
       AC_MSG_ERROR("no valid match for $2 found")
   else
       AC_MSG_RESULT("$INTEGER_SIZE_TYPE")
   fi
])

dnl-------------------------------------------------------------------------dnl
dnl DETERMINE C++ DATA TYPE FOR HOST_FLOAT
dnl
dnl DETERMINE C++ DATA TYPE FOR A GIVEN FLOAT SIZE
dnl eg. AC_DETERMINE_FLOAT(8) sets the variable FLOAT_SIZE_TYPE to double
dnl-------------------------------------------------------------------------dnl

AC_DEFUN(AC_DETERMINE_FLOAT, [dnl

   AC_MSG_CHECKING("C++ data type of float of size $1 bytes")

   # set the language to C++
   AC_LANG_CPLUSPLUS

   check_floats='false'
   
   # check to see if regular float does the job
   AC_TRY_RUN([
       int main()
       {
	   int p = 1;
	   if (sizeof(float) == $1)
	       p = 0;
	   return p;
       }], [FLOAT_SIZE_TYPE='float'], 
	   [check_floats='true'])
      
   # check double type if float did not pass
   if test "${check_floats}" = true ; then

       AC_TRY_RUN([
	   int main()
	   {   
	       int p = 1;
	       if (sizeof(double) == $1)
		   p = 0;
	       return p;
	   }], [check_floats='false';
	        FLOAT_SIZE_TYPE='double'], [])

   fi  
     
   # check long double type if double did not pass
   if test "${check_floats}" = true ; then

       AC_TRY_RUN([
	   int main()
	   {   
	       int p = 1;
	       if (sizeof(long double) == $1})
		   p = 0;
	       return p;
	   }], [check_floats='false';
	        FLOAT_SIZE_TYPE='long double'], [])

   fi
       
   # error message because we haven't found a valid type
   if test "${check_floats}" = true ; then
       AC_MSG_RESULT("no match found")
       AC_MSG_ERROR("no valid match for HOST_FLOATS found")
   else
       AC_MSG_RESULT("$FLOAT_SIZE_TYPE")
   fi
])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_local.m4
dnl-------------------------------------------------------------------------dnl
