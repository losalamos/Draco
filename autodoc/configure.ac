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

# draco version
AC_DRACO_PREREQ
AC_DRACO_ENV

dnl
dnl doxygen setup
dnl

# check for and as sign the path to doxygen
AC_PATH_PROG(DOXYGEN_PATH, doxygen, null)   

AC_MSG_CHECKING([for available Draco packages])
doxygen_examples=''
COMPLINKS=''
adp=''
for item in `ls -1 ${prefix}/html`; do
   if test -d ${prefix}/html/${item}; then
      dirname=`basename ${item}`
      adp="${adp} ${dirname}"
      COMPLINKS="${COMPLINKS} <li><a href=\"${dirname}/index.html\">${dirname}</a></li>"
   fi
done
AC_MSG_RESULT(${adp:-none})

# add config/doc to directories
doxygen_input=`pwd`
doxygen_input="${srcdir}/../config/doc ${srcdir} ${doxygen_input}"

# check if we should build latex sources
latex_yes_no='NO'
if test "${enable_latex_doc:=no}" = yes ; then
   latex_yes_no='YES'
fi

# dotpath (for draco levelization diagram)
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
AC_SUBST(COMPLINKS)

AC_CONFIG_FILES([Makefile doxygen_config header.html footer.html mainpage.dcc])
AC_OUTPUT

dnl-------------------------------------------------------------------------dnl
dnl			 end of draco/src/autodoc configure.in
dnl-------------------------------------------------------------------------dnl
