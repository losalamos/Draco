#!/bin/bash

#--------------------------------------------------------------------------------------------------#
# Provide the canonicalize filename (physical filename with out any symlinks) like the GNU version
# readlink with the -f option regardless of the version of readlink (GNU or BSD).
#
# This file is part of a set of unofficial pre-commit hooks available at github.
# Link:    https://github.com/githubbrowser/Pre-commit-hooks
# Contact: David Martin, david.martin.mailbox@googlemail.com
#
# Canonicalize by recursively following every symlink in every component of the specified filename.
# This should reproduce the results of the GNU version of readlink with the -f option.
#
# Reference: http://stackoverflow.com/questions/1055671/how-can-i-get-the-behavior-of-gnus-readlink-f-on-a-mac
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
# There should be no need to change anything below this line.
#--------------------------------------------------------------------------------------------------#
canonicalize_filename () {
  local target_file
  local physical_directory
  local result
  target_file=$1
  # $(echo $1 | head -n 1)
  physical_directory=""
  result=""

  # Need to restore the working directory after work.
  pushd "$(pwd)" &> /dev/null || \
    (echo "canonicalize_filename error 1: cannot access PWD" && return 1)

  cd "$(dirname "$target_file")" &> /dev/null || \
    (echo "canonicalize_filename error 2: cannot access target_file PWD" && return 2)
  target_file=$(basename "$target_file")

  # Iterate down a (possible) chain of symlinks
  while [ -L "$target_file" ]; do
    target_file=$(readlink "$target_file")
    cd "$(dirname "$target_file")" || \
      (echo "canonicalize_filename error 3: cannot access target_file PWD" && return 3)
    target_file=$(basename "$target_file")
  done

  # Compute the canonicalized name by finding the physical path for the directory we're in and
  # appending the target file.
  physical_directory=$(pwd -P)
  result="$physical_directory"/"$target_file"

  # restore the working directory after work.
  popd > /dev/null || (echo "canonicalize_filename error 4: cannot popd" && return 4)

  echo "${result}"
}

#--------------------------------------------------------------------------------------------------#
# End canonicalize_filename.sh
#--------------------------------------------------------------------------------------------------#
