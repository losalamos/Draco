#!/bin/sh

# File  : format-f90
# Author: Kelly Thompson, kgt@lanl.gov
# Date  : Wednesday, Sep 05, 2018, 16:01 pm

# Description:

# This small script will apply Emacs f90 formatting to the specified file.

# Usage:
#   ./format-f90.sh [-h] [file.f90]
#
# Args:
#   -h          Help message
#   -i          Don't do anything, just echo the commands
#-------------------------------------------------------------------------------

# Where am I?
export rscriptdir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" )
if ! [[ -d $rscriptdir ]]; then
  export rscriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
fi

# Defaults
EMACS=emacs

print_use()
{
  echo " "
  echo "This small script will apply Emacs formatting to the specified file."
  echo " "
  echo "Usage:"
  echo "  ./format-f90.sh [-hi] -f file.f90"
  echo " "
  echo "Args:"
  echo "  -f file.f90 The file you want to process."
  echo "  -h          Help message"
  echo "  -i          Only echo the commands, do not actually do anything."
}

# Parse arguments
f90source=""
while getopts "f:hi" opt; do
  case $opt in
    f) f90source=${OPTARG} ;;
    i) echoonly=echo ;;
    h) print_use; exit 0 ;;
    \?) echo "" ;echo "invalid option: -$OPTARG"; print_use; exit 1 ;;
    :)  echo "" ;echo "option -$OPTARG requires an argument."; print_use; exit 1 ;;

  esac
done

if ! [[ $f90source ]]; then
  echo -e "\nERROR: You must provide at least one file to format."
  print_use
  exit 1
fi

if [[ $echoonly ]]; then
  $echoonly $EMACS -batch ${f90source} -l ${rscriptdir}/format-f90.el -f emacs-format-f90-sources
else
  $EMACS -batch ${f90source} -l ${rscriptdir}/format-f90.el -f emacs-format-f90-sources &> /dev/null
fi
