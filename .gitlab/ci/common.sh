#!/bin/bash -l
# -*- Mode: sh -*-
#--------------------------------------------------------------------------------------------------#
# File  : .gitlab/ci/common.sh
# Date  : Tuesday, May 31, 2016, 14:48 pm
# Author: Kelly Thompson
# Note  : Copyright (C) 2016-2021, Triad National Security, LLC., All rights are reserved.
#--------------------------------------------------------------------------------------------------#
#
# Summary: Misc bash functions useful during development of code.

shopt -s extglob # Allow variable as case condition

#--------------------------------------------------------------------------------------------------#
# Helpful functions
#--------------------------------------------------------------------------------------------------#

function dracohelp()
{
echo -e "Bash functions defined by Draco:\n\n"
echo -e "\nUse 'type <function>' to print the full content of any function.\n"
echo -e "\ndie        - Print an error message and exit."
echo -e "\nrun        - echo a command and then evaluate it."
echo -e "\nfn_exists  - return true if the named bash function is defined."
echo -e "\nversion_gt - compare versions"
echo -e "\n"
}

#--------------------------------------------------------------------------------------------------#
# Function Definitions
#--------------------------------------------------------------------------------------------------#

# Print an error message and exit.  e.g.: cd $dir || die "can't change dir to $dir".
function die () { echo " "; echo "FATAL ERROR: $1"; exit 1;}

# Echo a command and then run it.
function run ()
{
  echo "==> $1"; if [[ "${dry_run:-no}" == "no" ]]; then eval "$1"; fi
}

# Return 0 if provided name is a bash function.
function fn_exists ()
{
  type "$1" 2>/dev/null | grep -c 'is a function'
}

#--------------------------------------------------------------------------------------------------#
# Compare versions
#
# Example:
#
# EMACSVER=`emacs --version | head -n 1 | sed -e 's/.*Emacs //'`
# if `version_gt "24.0.0" $EMACSVER` ; then echo yes; fi
#--------------------------------------------------------------------------------------------------#
function version_gt()
{
  test "$(printf '%s\n' "$@" | sort -V | head -n 1)" != "$1";
}

#--------------------------------------------------------------------------------------------------#
# check whether the given file matches any of the set extensions
#
# Example:
#
# FILE_EXTS=".f90 .F90"
# FILE_ENDINGS_INCLUDE="CMakeLists.txt"
# FILE_ENDINGS_EXCLUDE="_f.h _f77.h _f90.h"
# for file in $modified_files; do
#   if ! matches_extension $file; then continue; fi
#   <do stuff>
# done
#--------------------------------------------------------------------------------------------------#
matches_extension() {
  local filename
  local extension
  local end
  local ext
  filename=$(basename "$1")
  extension=".${filename##*.}"
  for end in $FILE_ENDINGS_EXCLUDE; do [[ "$filename" == *"$end" ]] && return 1; done
  for end in $FILE_ENDINGS_INCLUDE; do [[ "$filename" == *"$end" ]] && return 0; done
  for ext in $FILE_EXTS; do [[ "$ext" == "$extension" ]] && return 0; done
  return 1
}

#--------------------------------------------------------------------------------------------------#
# End common.sh
#--------------------------------------------------------------------------------------------------#
