#!/bin/bash

# Run clang-format in the current directory and list locally modified
# files that are not compliant with the current coding standard (see
# .clang_format in the top level source directory.)

##---------------------------------------------------------------------------##
## Environment
##---------------------------------------------------------------------------##

# Enable job control
set -m

##---------------------------------------------------------------------------##
## Support functions
##---------------------------------------------------------------------------##
print_use()
{
    echo " "
    echo "Usage: ${0##*/} -d -n -t"
    echo " "
    echo "All arguments are optional."
    echo "  -d Diff mode only. Do not modifty files."
    echo "  -n Alias for -d."
    echo "  -t Run as a pre-commit check, print list of non-conformant files and return"
    echo "     with exit code = 1."
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
gcf=`which git-clang-format${cfver}`
cf=`which clang-format${cfver}`
if [[ ! ${gcf} ]]; then
   echo "ERROR: git-clang-format${cfver} was not found in your PATH."
   echo "pwd="
   pwd
   echo "which git-clang-format${cfver}"
   echo $gcf
   echo "which clang-format${cfver}"
   echo $cf
   echo "which git"
   which git
   exit 1
else
  echo "Using $gcf --binary $cf"
fi

ver=`${cf} --version`
echo " "
echo "--------------------------------------------------------------------------------"
echo "Checking modified code for style conformance..."
echo "  - using clang-format version $ver"
echo "  - using settings from this project's .clang_format configuration file."
echo " "

##---------------------------------------------------------------------------##
## Default values
##---------------------------------------------------------------------------##
pct_mode=0
diff_mode=0
allok=0

##---------------------------------------------------------------------------##
## Command options
##---------------------------------------------------------------------------##

while getopts ":dhnt" opt; do
case $opt in
d) diff_mode=1 ;;
h) print_use; exit 0 ;;
n) diff_mode=1 ;;
t) pct_mode=1 ;;
\?) echo "" ;echo "invalid option: -$OPTARG"; print_use; exit 1 ;;
:)  echo "" ;echo "option -$OPTARG requires an argument."; print_use; exit 1 ;;
esac
done

##---------------------------------------------------------------------------##
## Check mode
##---------------------------------------------------------------------------##

if test "${pct_mode}" = "1"; then

  # don't actually modify the files (compare to branch 'develop')
  cmd='${gcf} --binary ${cf} -f --diff --extensions hh,cc develop'
  echo "Running..."
  echo "   ${gcf} --binary ${cf} -f --diff --extensions hh,cc develop"
  echo " "
  result=`eval $cmd`
  allok=`echo $result | grep -c "did not modify"`
  # 2nd chance (maybe there are no files to check)
  if test $allok = 0; then
    allok=`echo $result | grep -c "no modified files"`
  fi

  if test $allok = 1; then
    echo "PASS: Changes conform to this project's style requirements."
  else
    echo "FAIL: some files do not conform to this project's style requirements:"
    echo " "
    # rerun the command to capture color output.
    eval $cmd
    exit 1
  fi

##---------------------------------------------------------------------------##
## Fix mode
##   no options --> fix the files by running clang-format
##   -d | -n    --> print diff of required changes.
##---------------------------------------------------------------------------##

else

  if test ${diff_mode} = 1; then
    cmd='${gcf} --binary ${cf} -f --diff --extensions hh,cc develop'
    result=`eval $cmd`
    echo "The following non-conformances were discovered. Rerun without -d|-n to"
    echo "automatically apply these changes:"
    echo " "
    # rerun command to capture color output.
    eval $cmd
  else
    result=`${gcf} --binary ${cf} -f --extensions hh,cc develop`
    nonconformantfilesfound=`echo $result | grep -c "changed files"`
    echo "The following files in your working directory were modified to meet the"
    echo "this project's style requirement:"
    echo " "
    echo $result
  fi

fi

##---------------------------------------------------------------------------##
## End check_style.sh
##---------------------------------------------------------------------------##
