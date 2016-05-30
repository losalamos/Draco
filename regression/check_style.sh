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
    echo "Usage: `basename $0` -d -n -t"
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
if ! test "${CLANG_FORMAT_VER}x" = "x"; then
  cfver="-${CLANG_FORMAT_VER}"
else
  cfver=""
fi
gcf=`which git-clang-format${cfver}`
if test "${gcf}notset" = "notset"; then
   echo "ERROR: git-clang-format${cfver} was not found in your PATH."
   echo "pwd="
   pwd
   echo "which clang-format${cfver}"
   which clang-format${cfver}
   echo "which git"
   which git
   echo "find clang-format${cfver}"
   find . -name clang-format${cfver}
   exit 1
fi

ver=`clang-format${cfver} --version`
echo " "
echo "--------------------------------------------------------------------------------"
echo "Checking modified code for style conformance..."
echo "  - using clang-format version $ver"
echo "  - using settings from Draco's .clang_format configuration file."
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
  cmd='git-clang-format${cfver} -f --diff --extensions hh,cc develop'
  result=`eval $cmd`
  allok=`echo $result | grep "did not modify" | wc -l`
  # 2nd chance (maybe there are no files to check)
  if test $allok = 0; then
    allok=`echo $result | grep "no modified files" | wc -l`
  fi

  if test $allok = 1; then
    echo "PASS: Changes conform to draco style requirements."
  else
    echo "FAIL: some files do not conform to draco style requirements:"
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
    cmd='git-clang-format${cfver} -f --diff --extensions hh,cc develop'
    result=`eval $cmd`
    echo "The following non-conformances were discovered. Rerun without -d|-n to"
    echo "automatically apply these changes:"
    echo " "
    # rerun command to capture color output.
    eval $cmd
  else
    result=`git-clang-format${cfver} -f --extensions hh,cc develop`
    nonconformantfilesfound=`echo $result | grep "changed files" | wc -l`
    echo "The following files in your working directory were modified to meet the draco"
    echo "style requirement:"
    echo " "
    echo $result
  fi

fi

##---------------------------------------------------------------------------##
## End check_style.sh
##---------------------------------------------------------------------------##
