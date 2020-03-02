#!/bin/bash

# Run clang-format in the current directory and list locally modified
# files that are not compliant with the current coding standard (see
# .clang_format in the top level source directory.)

##---------------------------------------------------------------------------##
## Environment
##---------------------------------------------------------------------------##

# Enable job control
set -m

# protect temp files
umask 0077

# load some common bash functions
export rscriptdir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" )
if ! [[ -d $rscriptdir ]]; then
  export rscriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
fi
if [[ -f $rscriptdir/common.sh ]]; then
  source $rscriptdir/common.sh
else
  echo " "
  echo "FATAL ERROR: Unable to locate Draco's bash functions: "
  echo "   looking for .../tools/common.sh"
  echo "   searched rscriptdir = $rscriptdir"
  exit 1
fi

##---------------------------------------------------------------------------##
## Support functions
##---------------------------------------------------------------------------##
print_use()
{
    echo " "
    echo "Usage: ${0##*/} -f -t"
    echo " "
    echo "All arguments are optional."
    echo "  -f Show diff and fix files (when possible)."
    echo -n "  -t Run as a pre-commit check, print list of non-conformant "
    echo "files and return"
    echo "     with exit code = 1 (implies -d)."
    echo " "
}

##---------------------------------------------------------------------------##
## Sanity Checks
##---------------------------------------------------------------------------##

# clang-format must be in the PATH
if [[ ${CLANG_FORMAT_VER} ]]; then
  cfver="-${CLANG_FORMAT_VER}"
else
  cfver=""
fi
# Assume applications have version postfix.
gcf=`which git-clang-format${cfver}`
cf=`which clang-format${cfver}`
# if not found, try to find applications w/o version postfix.
if ! [[ -f ${gcf} ]]; then
  gcf=`which git-clang-format`
fi
if ! [[ -f ${cf} ]]; then
  gcf=`which clang-format`
fi
# if still not found, abort.
if [[ ! ${gcf} ]]; then
   echo "ERROR: git-clang-format${cfver} was not found in your PATH."
   echo "pwd="
   pwd
   echo "which git-clang-format${cfver}"
   echo $gcf
   exit 1
#else
#  echo "Using $gcf --binary $cf"
fi
if [[ ! ${cf} ]]; then
   echo "ERROR: clang-format${cfver} was not found in your PATH."
   echo "pwd="
   pwd
   echo "which clang-format${cfver}"
   echo $cf
   echo "which git"
   which git
   exit 1
fi

##---------------------------------------------------------------------------##
## Default values
##---------------------------------------------------------------------------##
pct_mode=0     # pre-commit test mode.
fix_mode=0     # show diffs AND modify code
foundissues=0  # 0 == ok

##---------------------------------------------------------------------------##
## Command options
##---------------------------------------------------------------------------##

while getopts ":fht" opt; do
case $opt in
f) fix_mode=1 ;; # also modify code (as possible)
h) print_use; exit 0 ;;
t) pct_mode=1    # pre-commit test
   fix_mode=0
   ;;
\?) echo "" ;echo "invalid option: -$OPTARG"; print_use; exit 1 ;;
:)  echo "" ;echo "option -$OPTARG requires an argument."; print_use; exit 1 ;;
esac
done

##---------------------------------------------------------------------------##
## Check mode (Test C++ code with git-clang-format)
##---------------------------------------------------------------------------##

ver=`${cf} --version`
echo -ne "\n------------------------------------------------------------"
echo "--------------------"
echo "Checking modified C/C++ code for style conformance..."
#echo "  - using clang-format version $ver"
#echo "  - using settings from this project's .clang_format configuration file."
echo " "

patchfile_c=$(mktemp /tmp/gcf.patch.XXXXXXXX)

# don't actually modify the files (originally we compared to branch 'develop',
# but let's try ORIG_HEAD or maybe use CI variables like TRAVIS_BRANCH or
# CI_MERGE_REQUEST_TARGET_BRANCH_NAME).
run "git branch -a"
echo "TRAVIS_BRANCH = $TRAVIS_BRANCH"
echo "CI_MERGE_REQUEST_TARGET_BRANCH_NAME = $CI_MERGE_REQUEST_TARGET_BRANCH_NAME"
target_branch=develop
if [[ -n ${TRAVIS_BRANCH} ]]; then
  target_branch=${TRAVIS_BRANCH}
elif [[ -n ${CI_MERGE_REQUEST_TARGET_BRANCH_NAME} ]]; then
  target_branch=${CI_MERGE_REQUEST_TARGET_BRANCH_NAME}
fi
echo "Looking at code changes compared to target branch = $target_branch"
cmd="${gcf} --binary ${cf} -f --diff --extensions hh,cc,cu $target_branch"
run "${cmd}" &> $patchfile_c

# if the patch file has the string "no modified files to format", the check
# passes.
if [[ `grep -c "no modified files" ${patchfile_c}` != 0 ]] || \
   [[ `grep -c "clang-format did not modify any files" ${patchfile_c}` != 0 ]]; then
  echo -n "PASS: Changes to C++ sources conform to this project's style "
  echo "requirements."
else
  foundissues=1
  echo -n "FAIL: some C++ files do not conform to this project's style "
  echo "requirements:"
  # Modify files, if requested
  if [[ ${fix_mode} == 1 ]]; then
    echo -e "      The following patch has been applied to your file.\n"
    run "git apply $patchfile_c"
    cat $patchfile_c
  else
    echo -ne "      run ${0##*/} with option -f to automatically apply this "
    echo -e "patch.\n"
    cat $patchfile_c
  fi
fi
rm -f $patchfile_c

##---------------------------------------------------------------------------##
## Check mode (Test F90 code indentation with emacs and bash)
##---------------------------------------------------------------------------##

# Defaults ----------------------------------------
EMACS=`which emacs`
if [[ $EMACS ]]; then
  EMACSVER=`$EMACS --version | head -n 1 | sed -e 's/.*Emacs //'`
  if `version_gt "24.0.0" $EMACSVER` ; then
    echo "WARNING: Your version of emacs is too old. Expecting v 24.0+."
    echo "         pre-commit-hook partially disabled (f90 indentation)"
    unset EMACS
  fi
fi

# staged files
modifiedfiles=`git diff --name-only --cached`
# unstaged files
modifiedfiles="${modifiedfiles} `git diff --name-only`"
# all files
modifiedfiles=`echo ${modifiedfiles} | sort -u`

# file types to parse. Only effective when PARSE_EXTS is true.
FILE_EXTS=".f90 .F90"

# file endings for files to exclude from parsing when PARSE_EXTS is true.
FILE_ENDINGS="_f.h _f77.h _f90.h"

if [[ -x $EMACS ]]; then

  echo -ne "\n----------------------------------------------------------------"
  echo "----------------"
  echo -e "Checking modified F90 code for style conformance (indentation)..\n"

  patchfile_f90=$(mktemp /tmp/emf90.patch.XXXXXXXX)

  # Loop over all modified F90 files.  Create one patch containing all changes
  # to these files
  for file in $modifiedfiles; do

    # ignore file if we do check for file extensions and the file does not match
    # any of the extensions specified in $FILE_EXTS
    if ! matches_extension "$file"; then continue; fi

    file_nameonly=`basename $file`
    tmpfile1=/tmp/f90-format-$file_nameonly
    cp -f $file $tmpfile1
    $EMACS -batch ${tmpfile1} -l ${rscriptdir}/../environment/git/f90-format.el \
      -f emacs-format-f90-sources &> /dev/null
    # color output is possible if diff -version >= 3.4 with option `--color`
    diff -u "${file}" "${tmpfile1}" | \
      sed -e "1s|--- |--- a/|" -e "2s|+++ ${tmpfile1}|+++ b/${file}|" \
          >> "$patchfile_f90"
    rm $tmpfile1

  done

  # If the patch file is size 0, then no changes are needed.
  if [[ -s "$patchfile_f90" ]]; then
    foundissues=1
    echo -n "FAIL: some F90 files do not conform to this project's style "
    echo "requirements:"
    # Modify files, if requested
    if [[ ${fix_mode} == 1 ]]; then
      echo -e "      The following patch has been applied to your file.\n"
      run "git apply $patchfile_f90"
      cat $patchfile_f90
    else
      echo -ne "      run ${0##*/} with option -f to automatically apply this "
      ehco -e "patch.\n"
      cat $patchfile_f90
    fi
  else
    echo -n "PASS: Changes to F90 sources conform to this project's style "
    echo "requirements."
  fi
  rm -f $patchfile_f90

fi

##---------------------------------------------------------------------------##
## Check mode (Test F90 code line length with bash)
##---------------------------------------------------------------------------##

echo -ne "\n----------------------------------------------------------------"
echo "----------------"
echo -e "Checking modified F90 code for style conformance (line length)..\n"

tmpfile2=$(mktemp /tmp/f90-format-line-len.XXXXXXXX)
# Loop over all modified F90 files.  Create one patch containing all changes
# to these files
for file in $modifiedfiles; do

  # ignore file if we do check for file extensions and the file does not match
  # any of the extensions specified in $FILE_EXTS
  if ! matches_extension "$file"; then continue; fi

  header_printed=0
  lindeno=0

  cat "$file" | while read line; do
    let lineno++
    # Exceptions:
    # - Long URLs
    exception=`echo $line | grep -i -c http`
    if [[ ${#line} -gt 80 && ${exception} == 0 ]]; then
      if [[ ${header_printed} == 0 ]]; then
        echo -e "\nFile: ${file} [code line too long]\n" >> $tmpfile2
        echo "  line   length content" >> $tmpfile2
        echo -n "  ------ ------ -------------------------------" >> $tmpfile2
        echo "-------------------------------------------------" >> $tmpfile2
        header_printed=1
      fi
      printf "  %-6s %-6s %s\n" "${lineno}" "${#line}" "${line}" >> $tmpfile2
    fi
    # reset exception flag
    exception=0
  done
done

# If there are issues, report them
if [[ `cat $tmpfile2 | wc -l` -gt 0 ]]; then
    foundissues=1
    echo -ne "FAIL: some F90 files do not conform to this project's style "
    echo "requirements:\n"
    cat $tmpfile2
    echo -ne "\nPlease reformat lines listed above to fit into 80 columns and "
    echo -ne "attempt running\n${0##*/} again. These issues cannot be fixed "
    echo "with the -f option."
fi

#------------------------------------------------------------------------------#
# Done
#------------------------------------------------------------------------------#

# Return code: 0==ok,1==bad
exit $foundissues

##---------------------------------------------------------------------------##
## End check_style.sh
##---------------------------------------------------------------------------##
