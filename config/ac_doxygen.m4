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
dnl AC_SET_DEFAULT_OUTPUT
dnl-------------------------------------------------------------------------dnl
#
# Set the default location for doxygen output
#
AC_DEFUN([AC_SET_DEFAULT_OUTPUT], [dnl
   if test ${doxygen_output_top} = DEFAULT; then
       AC_SUBST(doxygen_output_top, "${prefix}/documentation")
   fi
])


dnl-------------------------------------------------------------------------dnl
dnl AC_AUTODOC_COMPONENT_TAGS
dnl
dnl   Collect tagfiles for within-package component dependencies
dnl-------------------------------------------------------------------------dnl
#
# Build a list of tagfiles for other components of the same package
# and the _relative_ locations of the autodoc directories that they
# refer to.
#
# The relative path between documentation for components in the
# same package is "../component"
#
# These components are specified in AC_NEEDS_LIBS, and are stored
# in variable DEPENDENT_COMPONENTS. 
#
AC_DEFUN([AC_AUTODOC_COMPONENT_TAGS], [dnl

   components=''
   AC_MSG_CHECKING([for Doxygen component dependencies])
   for comp in ${DEPENDENT_COMPONENTS}; do
       components="${components} ${comp}"
       TAGFILES="${TAGFILES} ${doxygen_output_top}/${comp}.tag"
       DOXYGEN_TAGFILES="${DOXYGEN_TAGFILES} \"${doxygen_output_top}/${comp}.tag = ../${comp}\""
   done
   AC_MSG_RESULT([${components}])

   AC_SUBST(TAGFILES)
   AC_SUBST(DOXYGEN_TAGFILES)

])

dnl-------------------------------------------------------------------------dnl
dnl AC_DRACO_AUTODOC
dnl
dnl  setup doxygen autodoc directories for COMPONENTS within a package
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_DRACO_AUTODOC], [dnl

   # Define some package-level directories
   header_dir=${package_top_srcdir}/autodoc/html
   autodoc_dir=${doxygen_input}/autodoc

   # For a component, the doxygen input is the srcdir and the examples
   # are in the tests
   doxygen_input=`cd ${srcdir}; pwd`
   doxygen_examples=${doxygen_input}/test

   # Get the default output location
   AC_SET_DEFAULT_OUTPUT

   # Set the package-level html output location
   package_html=${doxygen_output_top}/html

   # The local dir is different from the current dir.
   localdir=`pwd`/autodoc

   # Set the component output locations.
   doxygen_html_output="${doxygen_output_top}/html/${package}"
   doxygen_latex_output="${doxygen_output_top}/latex/${package}"

   # compute relative paths from localdir
   adl_COMPUTE_RELATIVE_PATHS([\
      localdir:doxygen_examples:rel_doxygen_examples \
      localdir:doxygen_input:rel_doxygen_input\
      doxygen_html_output:package_html:rel_package_html\
   ])

   # add autodoc directory under the source to doxygen input
   rel_doxygen_input="$rel_doxygen_input $rel_doxygen_input/autodoc"

   # Get tags for other components in this package which this
   # component depends on
   AC_AUTODOC_COMPONENT_TAGS

   # XXX We will need to expand this to handle tag files in other
   # packages too.

   # find the release number
   number=$1
   AC_MSG_CHECKING("component release number")
   AC_MSG_RESULT($number)

   AC_SUBST(doxygen_input)
   AC_SUBST(doxygen_examples)
   AC_SUBST(doxygen_output_top)

   AC_SUBST(doxygen_html_output)
   AC_SUBST(doxygen_latex_output)

   AC_SUBST(rel_doxygen_examples)
   AC_SUBST(rel_doxygen_input)
   AC_SUBST(rel_package_html)

   AC_SUBST(dotpath)
   AC_SUBST(number)

   AC_SUBST(header_dir)
   AC_SUBST(autodoc_dir)

   AC_CONFIG_FILES([autodoc/Makefile:../../config/Makefile.autodoc.in \
                    autodoc/doxygen_config:../../config/doxygen_config.in \
                    autodoc/header.html:../../autodoc/html/header.html.in \
                    autodoc/footer.html:../../autodoc/html/footer.html.in ])

])

dnl-------------------------------------------------------------------------dnl
dnl AC_PACKAGE_AUTODOC
dnl
dnl  setup doxygen autodoc directories for a PACKAGE
dnl-------------------------------------------------------------------------dnl

AC_DEFUN([AC_PACKAGE_AUTODOC], [dnl

   AC_SET_DEFAULT_OUTPUT

   #
   # For package-level documentation, the only doxygen sources are in
   # the build tree, aka, the current directory.
   #
   #   doxygen_input=$builddir
   #   doxygen_input="${srcdir}/../config/doc ${srcdir} ${doxygen_input}"

   doxygen_input=`pwd`

   localdir=`pwd`

   #
   # compute relative paths from localdir
   #
   adl_COMPUTE_RELATIVE_PATHS([\
      localdir:doxygen_input:rel_doxygen_input \
   ])


   doxygen_html_output="${doxygen_output_top}/html/"
   doxygen_latex_output="${doxygen_output_top}/latex/"


   #
   # XXX Need to change COMPLINKS to generic doxygen list instead of
   # HTML for Latex compatability. Let doxygen insert the links
   #
   AC_MSG_CHECKING([for documented sub-components of this package])
   COMP_LINKS=''
   components=''
   for item in `ls -1 ${package_top_srcdir}/src`; do
      if test -d ${package_top_srcdir}/src/${item}/autodoc; then
         dirname=`basename ${item}`
         components="${components} ${dirname}"
         COMP_LINKS="${COMP_LINKS} <li><a href=\"${dirname}/index.html\">${dirname}</a></li>"
      fi
   done
   AC_MSG_RESULT(${components:-none})
   COMP_LINKS="<ul> $COMP_LINKS </ul>"

   #
   # XXX TO DO: Add links to dependent packages on this page.
   #
   PACKAGE_LINKS="<ul> </ul>"


   # XXX Need to build a list of tags for the doxygenated
   # components of the package.
   TAGFILES=''
   DOXYGEN_TAGFILES=''
   AC_SUBST(TAGFILES)
   AC_SUBST(DOXYGEN_TAGFILES)

   header_dir=${srcdir}/html
   AC_SUBST(header_dir)

   autodoc_dir=${srcdir}
   AC_SUBST(autodoc_dir)

   AC_SUBST(rel_doxygen_input)
   AC_SUBST(rel_doxygen_examples, '')
   AC_SUBST(doxygen_html_output)
   AC_SUBST(doxygen_latex_output)

   AC_SUBST(PACKAGE_LINKS)
   AC_SUBST(COMP_LINKS)

   AC_CONFIG_FILES([doxygen_config:../config/doxygen_config.in])
   AC_CONFIG_FILES([Makefile:../config/Makefile.autodoc.in])
   AC_CONFIG_FILES([header.html:html/header.html.in])
   AC_CONFIG_FILES([footer.html:html/footer.html.in])
   AC_CONFIG_FILES([mainpage.dcc])

])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_doxygen.m4
dnl-------------------------------------------------------------------------dnl

