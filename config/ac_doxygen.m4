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

   adl_COMPUTE_RELATIVE_PATH([localdir],[doxygen_output],[rel_doxygen_output])
   adl_COMPUTE_RELATIVE_PATH([localdir],[doxygen_examples],[rel_doxygen_examples])
   adl_COMPUTE_RELATIVE_PATH([localdir],[doxygen_input],[rel_doxygen_input])

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

dnl http://www.gnu.org/software/ac-archive/htmldoc/relpaths.html
dnl

AC_DEFUN([adl_COMPUTE_RELATIVE_PATHS],
[for _lcl_i in $1; do
  _lcl_from=\[$]`echo "[$]_lcl_i" | sed 's,:.*$,,'`
  _lcl_to=\[$]`echo "[$]_lcl_i" | sed 's,^[[^:]]*:,,' | sed 's,:[[^:]]*$,,'`
  _lcl_result_var=`echo "[$]_lcl_i" | sed 's,^.*:,,'`
  adl_RECURSIVE_EVAL([[$]_lcl_from], [_lcl_from])
  adl_RECURSIVE_EVAL([[$]_lcl_to], [_lcl_to])
  _lcl_notation="$_lcl_from$_lcl_to"
  adl_NORMALIZE_PATH([_lcl_from],['/'])
  adl_NORMALIZE_PATH([_lcl_to],['/'])
  adl_COMPUTE_RELATIVE_PATH([_lcl_from], [_lcl_to], [_lcl_result_tmp])
  adl_NORMALIZE_PATH([_lcl_result_tmp],["[$]_lcl_notation"])
  eval $_lcl_result_var='[$]_lcl_result_tmp'
done])

## Note:
## *****
## The following helper macros are too fragile to be used out
## of adl_COMPUTE_RELATIVE_PATHS (mainly because they assume that
## paths are normalized), that's why I'm keeping them in the same file.
## Still, some of them maybe worth to reuse.

dnl adl_COMPUTE_RELATIVE_PATH(FROM, TO, RESULT)
dnl ===========================================
dnl Compute the relative path to go from $FROM to $TO and set the value
dnl of $RESULT to that value.  This function work on raw filenames
dnl (for instead it will considerate /usr//local and /usr/local as
dnl two distinct paths), you should really use adl_COMPUTE_REALTIVE_PATHS
dnl instead to have the paths sanitized automatically.
dnl
dnl For instance:
dnl    first_dir=/somewhere/on/my/disk/bin
dnl    second_dir=/somewhere/on/another/disk/share
dnl    adl_COMPUTE_RELATIVE_PATH(first_dir, second_dir, first_to_second)
dnl will set $first_to_second to '../../../another/disk/share'.
AC_DEFUN([adl_COMPUTE_RELATIVE_PATH],
[adl_COMPUTE_COMMON_PATH([$1], [$2], [_lcl_common_prefix])
adl_COMPUTE_BACK_PATH([$1], [_lcl_common_prefix], [_lcl_first_rel])
adl_COMPUTE_SUFFIX_PATH([$2], [_lcl_common_prefix], [_lcl_second_suffix])
$3="[$]_lcl_first_rel[$]_lcl_second_suffix"])

dnl adl_COMPUTE_COMMON_PATH(LEFT, RIGHT, RESULT)
dnl ============================================
dnl Compute the common path to $LEFT and $RIGHT and set the result to $RESULT.
dnl
dnl For instance:
dnl    first_path=/somewhere/on/my/disk/bin
dnl    second_path=/somewhere/on/another/disk/share
dnl    adl_COMPUTE_COMMON_PATH(first_path, second_path, common_path)
dnl will set $common_path to '/somewhere/on'.
AC_DEFUN([adl_COMPUTE_COMMON_PATH],
[$3=''
_lcl_second_prefix_match=''
while test "[$]_lcl_second_prefix_match" != 0; do
  _lcl_first_prefix=`expr "x[$]$1" : "x\([$]$3/*[[^/]]*\)"`
  _lcl_second_prefix_match=`expr "x[$]$2" : "x[$]_lcl_first_prefix"`
  if test "[$]_lcl_second_prefix_match" != 0; then
    if test "[$]_lcl_first_prefix" != "[$]$3"; then
      $3="[$]_lcl_first_prefix"
    else
      _lcl_second_prefix_match=0
    fi
  fi
done])

dnl adl_COMPUTE_SUFFIX_PATH(PATH, SUBPATH, RESULT)
dnl ==============================================
dnl Substrack $SUBPATH from $PATH, and set the resulting suffix
dnl (or the empty string if $SUBPATH is not a subpath of $PATH)
dnl to $RESULT.
dnl
dnl For instace:
dnl    first_path=/somewhere/on/my/disk/bin
dnl    second_path=/somewhere/on
dnl    adl_COMPUTE_SUFFIX_PATH(first_path, second_path, common_path)
dnl will set $common_path to '/my/disk/bin'.
AC_DEFUN([adl_COMPUTE_SUFFIX_PATH],
[$3=`expr "x[$]$1" : "x[$]$2/*\(.*\)"`])

dnl adl_COMPUTE_BACK_PATH(PATH, SUBPATH, RESULT)
dnl ============================================
dnl Compute the relative path to go from $PATH to $SUBPATH, knowing that
dnl $SUBPATH is a subpath of $PATH (any other words, only repeated '../'
dnl should be needed to move from $PATH to $SUBPATH) and set the value
dnl of $RESULT to that value.  If $SUBPATH is not a subpath of PATH,
dnl set $RESULT to the empty string.
dnl
dnl For instance:
dnl    first_path=/somewhere/on/my/disk/bin
dnl    second_path=/somewhere/on
dnl    adl_COMPUTE_BACK_PATH(first_path, second_path, back_path)
dnl will set $back_path to '../../../'.
AC_DEFUN([adl_COMPUTE_BACK_PATH],
[adl_COMPUTE_SUFFIX_PATH([$1], [$2], [_lcl_first_suffix])
$3=''
_lcl_tmp='xxx'
while test "[$]_lcl_tmp" != ''; do
  _lcl_tmp=`expr "x[$]_lcl_first_suffix" : "x[[^/]]*/*\(.*\)"`
  if test "[$]_lcl_first_suffix" != ''; then
     _lcl_first_suffix="[$]_lcl_tmp"
     $3="../[$]$3"
  fi
done])

dnl adl_RECURSIVE_EVAL(VALUE, RESULT)
dnl =================================
dnl Interpolate the VALUE in loop until it doesn't change,
dnl and set the result to $RESULT.
dnl WARNING: It's easy to get an infinite loop with some unsane input.
AC_DEFUN([adl_RECURSIVE_EVAL],
[_lcl_receval="$1"
$2=`(test "x$prefix" = xNONE && prefix="$ac_default_prefix"
     test "x$exec_prefix" = xNONE && exec_prefix="${prefix}"
     _lcl_receval_old=''
     while test "[$]_lcl_receval_old" != "[$]_lcl_receval"; do
       _lcl_receval_old="[$]_lcl_receval"
       eval _lcl_receval="\"[$]_lcl_receval\""
     done
     echo "[$]_lcl_receval")`])

dnl-------------------------------------------------------------------------dnl
dnl end of ac_doxygen.m4
dnl-------------------------------------------------------------------------dnl

