dnl-------------------------------------------------------------------------dnl
dnl draco/src/autodoc configure.in
dnl Thomas M. Evans
dnl Fri Sep  3 13:04:05 1999
dnl $Id$
dnl-------------------------------------------------------------------------dnl
dnl @> configure.in for autodoc
dnl-------------------------------------------------------------------------dnl

AC_INIT
AC_CONFIG_SRCDIR([doxygen_config.in])
AC_CONFIG_AUX_DIR(../config)

dnl
dnl doxygen setup
dnl

# check for and as sign the path to doxygen
AC_PATH_PROG(DOXYGEN_PATH, doxygen, null)   

# path to package directories which are sources for doxygen
doxygen_src=${srcdir}/../src

# find sources
AC_MSG_CHECKING("directories to doxygenize")

dirs=''
if test -d ../src ; then
   dirs=`ls ../src`
fi

doxygen_input=''
doxygen_examples=''
doxygen_dirs=''

for dir in $dirs;
do
   if test -d ../src/$dir ; then
       doxygen_dirs="${dir} ${doxygen_dirs}"
       doxygen_input="${doxygen_src}/${dir} ${doxygen_input}"
       doxygen_examples="${doxygen_src}/${dir}/test ${doxygen_examples}"
   fi
done

# add config/doc to directories
doxygen_dirs="config ${doxygen_dirs}"
doxygen_input="${srcdir}/../config/doc ${doxygen_input}"

# add autodoc to directories
doxygen_dirs="autodoc ${doxygen_dirs}"
doxygen_input="${srcdir} ${doxygen_input}"

AC_MSG_RESULT($doxygen_dirs)

# check if we should build latex sources
latex_yes_no='NO'
if test "${enable_latex_doc:=no}" = yes ; then
   latex_yes_no='YES'
fi

# dotpath (for draco levelization diagram
dotpath="${srcdir}/../doc/design"

# find the release number
number='null'
AC_MSG_CHECKING("draco release number")
number=`/bin/sh ../tools/draco --version`
AC_MSG_RESULT($number)

AC_SUBST(doxygen_input)
AC_SUBST(doxygen_examples)
AC_SUBST(latex_yes_no)
AC_SUBST(dotpath)
AC_SUBST(number)

AC_CONFIG_FILES([Makefile doxygen_config header.html footer.html])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/autodoc configure.in
dnl-------------------------------------------------------------------------dnl
