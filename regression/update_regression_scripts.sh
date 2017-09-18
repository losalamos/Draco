#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/update_regression_scripts.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# switch to group 'ccsrad' and set umask
if [[ $(id -gn) != ccsrad ]]; then
  exec sg ccsrad "$0 $*"
fi
umask 0002

# Locate the directory that this script is located in:
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Redirect all output to a log file.
timestamp=`date +%Y%m%d-%H%M`
target="`uname -n | sed -e s/[.].*//`"
logdir="$( cd $scriptdir/../../logs && pwd )"
logfile=$logdir/update_regression_scripts-$target-$timestamp.log
exec > $logfile
exec 2>&1

# import some bash functions
source $scriptdir/scripts/common.sh

echo -e "Executing $0 $*...\n"
echo "Group: `id -gn`"
echo -e "umask: `umask` \n"

#------------------------------------------------------------------------------#

# Per machine setup
case ${target} in
  darwin-fe* | cn[0-9]*)
    REGDIR=/usr/projects/draco/regress
    keychain=keychain-2.7.1
    VENDOR_DIR=/usr/projects/draco/vendors
    # personal copy of ssh-agent.
    export PATH=$HOME/bin:$PATH
    export http_proxy=http://proxyout.lanl.gov:8080;
    export https_proxy=$http_proxy;
    export HTTP_PROXY=$http_proxy;
    export HTTPS_PROXY=$http_proxy;
    export no_proxy="localhost,127.0.0.1,.lanl.gov";
    export NO_PROXY=$no_proxy;
    ;;
  ccscs*)
    REGDIR=/scratch/regress
    VENDOR_DIR=/scratch/vendors
    keychain=keychain-2.8.2
    ;;
  ml-*)
    REGDIR=/usr/projects/jayenne/regress
    keychain=keychain-2.7.1
    VENDOR_DIR=/usr/projects/draco/vendors
    ;;
  *)
    REGDIR=/scratch/regress
    ;;
esac

# Load some identities used for accessing gitlab.
if [[ -f $HOME/.ssh/id_rsa ]]; then
  MYHOSTNAME="`uname -n`"
  run "$VENDOR_DIR/$keychain/keychain $HOME/.ssh/id_rsa"
  if [[ -f $HOME/.keychain/$MYHOSTNAME-sh ]]; then
    run "source $HOME/.keychain/$MYHOSTNAME-sh"
  else
    echo "Error: could not find $HOME/.keychain/$MYHOSTNAME-sh"
  fi
fi

# ---------------------------------------------------------------------------- #
# Update the regression script directories
# ---------------------------------------------------------------------------- #

# Deal with proxy stuff on darwin
case ${target} in
  darwin-fe* | cn[0-9]*)
    unset http_proxy;
    unset https_proxy;
    unset HTTP_PROXY;
    unset HTTPS_PROXY;
    unset no_proxy;
    unset NO_PROXY;
  ;;
esac

# Setup
if ! [[ -d $REGDIR ]]; then
  run "mkdir -p ${REGDIR}"
  run "chmod g+rwX,o-rwX ${REGDIR}"
  run "chmod g+s ${REGDIR}"
fi

# Draco/Jayenne/Capsaicin
projects="draco jayenne capsaicin"
for p in $projects; do
  echo -e "\nUpdating $REGDIR/$p..."
  if test -d ${REGDIR}/$p; then
    run "cd ${REGDIR}/$p; git pull"
  else
    case $p in
      draco)
        run "cd ${REGDIR}; git clone git@github.com:lanl/Draco.git $p"
        ;;
      *)
        run "cd ${REGDIR}; git clone git@gitlab.lanl.gov:${p}/${p}.git"
        ;;
    esac
  fi
done

#------------------------------------------------------------------------------#
# Cleanup old files and directories
#------------------------------------------------------------------------------#
if [[ -d $logdir ]]; then
  echo -e "\nCleaning up old log files."
  run "cd $logdir"
  run "find . -mtime +14 -type f"
  run "find . -mtime +14 -type f -delete"
fi
if [[ -d $REGDIR/cdash ]]; then
  echo -e "\nCleaning up old builds."
  run "cd $REGDIR/cdash"
  run "find . -maxdepth 3 -mtime +14 -name 'Experimental*-pr*' -type d"
  run "find . -maxdepth 3 -mtime +14 -name 'Experimental*-pr*' -type d -exec rm -rf {} \;"
fi

echo -e "\n--------------------------------------------------------------------------------"
echo "All done."
echo "--------------------------------------------------------------------------------"

##---------------------------------------------------------------------------##
## End update_regression_scripts.sh
##---------------------------------------------------------------------------##
