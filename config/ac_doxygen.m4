dnl-------------------------------------------------------------------------dnl
dnl ac_doxygen.m4
dnl
dnl Macros to help setup doxygen autodoc directories.
dnl
dnl Kelly Thompson
dnl 2004/03/30 16:41:22
dnl 1999/02/04 01:56:19
dnl-------------------------------------------------------------------------dnl
##---------------------------------------------------------------------------##
## $Id$
##---------------------------------------------------------------------------##

dnl-------------------------------------------------------------------------dnl
dnl AC_DRACO_AUTODOC
dnl
dnl  setup doxygen autodoc directories.
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DRACO_AUTODOC], [dnl

   #
   # paths of package directories which are sources for doxygen
   #

   doxygen_input=`cd ${srcdir}; pwd`
   doxygen_examples=`cd ${srcdir}/test; pwd`
   localdir=`pwd`/autodoc

   # XXX This will change with new destination definition.
   doxygen_output_top="${prefix}/html"
   doxygen_output="${doxygen_output_top}/${package}"

   #
   # compute relative paths from localdir
   #

   adl_COMPUTE_RELATIVE_PATHS([\
      localdir:doxygen_output:rel_doxygen_output \
      localdir:doxygen_examples:rel_doxygen_examples \
      localdir:doxygen_input:rel_doxygen_input\
   ])

   #
   # add autodoc directory to doxygen input
   # 
   rel_doxygen_input="$rel_doxygen_input $rel_doxygen_input/autodoc"

   #
   # use relative paths for tag files also
   #
   components=''
   AC_MSG_CHECKING([for Doxygen component dependencies])
   for comp in ${DRACO_COMPONENTS}; do
       components="${components} ${comp}"
       TAGFILES="${TAGFILES} ${doxygen_output_top}/${comp}.tag"
       DOXYGEN_TAGFILES="${DOXYGEN_TAGFILES} \"${doxygen_output_top}/${comp}.tag = ../${comp}\""
   done
   AC_MSG_RESULT([${components}])

   # check if we should build latex sources
   latex_yes_no='NO'
   if test "${enable_latex_doc:=no}" = yes ; then
      latex_yes_no='YES'
   fi

   # find the release number
   number=$1
   AC_MSG_CHECKING("component release number")
   AC_MSG_RESULT($number)

   AC_SUBST(doxygen_output)
   AC_SUBST(doxygen_examples)
   AC_SUBST(doxygen_input)

   AC_SUBST(rel_doxygen_output)
   AC_SUBST(rel_doxygen_examples)
   AC_SUBST(rel_doxygen_input)

   AC_SUBST(latex_yes_no)
   AC_SUBST(dotpath)
   AC_SUBST(number)
   AC_SUBST(TAGFILES)
   AC_SUBST(DOXYGEN_TAGFILES)

   AC_CONFIG_FILES([autodoc/Makefile:../../config/Makefile.autodoc.in \
                    autodoc/doxygen_config:../../config/doxygen_config.in])

])


dnl-------------------------------------------------------------------------dnl
dnl end of ac_doxygen.m4
dnl-------------------------------------------------------------------------dnl

