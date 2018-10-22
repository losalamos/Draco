#!/bin/bash

# Install git pre-commit hooks. The root directory of the target local git
# repository is expected as a parameter.

# Example:
# .../draco/environment/git/install-hooks.sh $HOME/ccs-codes/jayenne

# This file is part of a set of unofficial pre-commit hooks available
# at github.
# Link:    https://github.com/githubbrowser/Pre-commit-hooks
# Contact: David Martin, david.martin.mailbox@googlemail.com

###########################################################
# CONFIGURATION:
# select which pre-commit hooks are going to be installed
#HOOKS="pre-commit pre-commit-compile pre-commit-uncrustify"
HOOKS="pre-commit pre-commit-clang-format"
###########################################################
# There should be no need to change anything below this line.

export rscriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [[ $1 ]]; then
  if [[ -d $1/.git ]]; then
    export dotgitdir="$( cd "$1/.git" && pwd )"
  else
    echo "ERROR: Cannot find $1/.git. Aborting"
    exit 1
  fi
else
  export dotgitdir="$( cd "$rscriptdir/../../.git" && pwd )"
fi

. "$(dirname -- "$0")/canonicalize_filename.sh"

# exit on error
set -e

# Absolute path to this script, e.g. /home/user/bin/foo.sh
SCRIPT="$(canonicalize_filename "$0")"

# Absolute path this script is in, e.g. /home/user/bin/
SCRIPTPATH="$(dirname -- "$SCRIPT")"

# copy hooks to the directory specified as parameter
copy_hooks() {
  echo "Copying hooks to destination directory:"
  for hook in $HOOKS
  do
    echo "Copying $hook to $dotgitdir/hooks/."
    cp -i -- "$SCRIPTPATH/$hook" "$dotgitdir/hooks/." || true
  done

  echo ""
  echo "Checking if hooks are executable:"
  for hook in $HOOKS
  do
    if [ ! -x "$dotgitdir/hooks/$hook" ] ; then
      chmod +x $dotgitdir/hooks/$hook
    else
      echo "$hook OK."
    fi
  done
}


echo ""
echo "Git pre-commit hook installation."
echo ""


if [ -d "$dotgitdir/hooks" ] ; then
  # create hooks subfolder if it does not yet exist
  mkdir -p -- "$dotgitdir/hooks"

  echo "Copying prerequisites."
  cp -i -- "$SCRIPTPATH/canonicalize_filename.sh" "$dotgitdir/hooks" || true
  echo ""
  copy_hooks
  echo ""
  echo "Finished installation."
else
  echo "Error: $dotgitdir does not exist."
  echo "Are you sure $dotgitdir is the root directory of your local git repository?"
fi
