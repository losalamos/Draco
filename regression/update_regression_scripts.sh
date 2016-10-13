#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/update_regression_scripts.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

umask 0002

target="`uname -n | sed -e s/[.].*//`"

# Locate the directory that this script is located in:
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# import some bash functions
source $scriptdir/scripts/common.sh

# Per machine setup
case ${target} in
  darwin-fe* | cn[0-9]*)
    REGDIR=/usr/projects/draco/regress
    keychain=keychain-2.7.1
    VENDOR_DIR=/usr/projects/draco/vendors
    # personal copy of ssh-agent.
    export PATH=$HOME/bin:$PATH
    if test -d /projects/opt/centos7/subversion/1.9.2/bin; then
      export PATH=/projects/opt/centos7/subversion/1.9.2/bin:$PATH
    fi
    SVN=`which svn`
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
    if ! [[ $MODULESHOME ]]; then
       source /usr/share/Modules/init/bash
    fi
    run "module load user_contrib subversion"
    SVN=`which svn`
    ;;
  ml-*)
    REGDIR=/usr/projects/jayenne/regress
    keychain=keychain-2.7.1
    VENDOR_DIR=/usr/projects/draco/vendors
    SVN=/usr/projects/hpcsoft/toss2/common/subversion/1.9.1/bin/svn
    ;;
  *)
    # module load user_contrib subversion
    SVN=/scratch/vendors/subversion-1.9.3/bin/svn
    REGDIR=/scratch/regress
    ;;
esac

# Load some identities used for accessing gitlab.
MYHOSTNAME="`uname -n`"
$VENDOR_DIR/$keychain/keychain $HOME/.ssh/cmake_dsa
if test -f $HOME/.keychain/$MYHOSTNAME-sh; then
    run "source $HOME/.keychain/$MYHOSTNAME-sh"
else
    echo "Error: could not find $HOME/.keychain/$MYHOSTNAME-sh"
fi

# ---------------------------------------------------------------------------- #
# Update the regression script directories
# ---------------------------------------------------------------------------- #

# Draco
echo " "
echo "Updating $REGDIR/draco..."
if ! test -d $REGDIR; then
  run "mkdir -p ${REGDIR}"
fi
if test -d ${REGDIR}/draco; then
  run "cd ${REGDIR}/draco; git pull"
else
  run "cd ${REGDIR}; git clone https://github.com/losalamos/Draco.git draco"
fi
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

# Jayenne
echo " "
echo "Updating $REGDIR/jayenne..."
if test -d ${REGDIR}/jayenne; then
  run "cd ${REGDIR}/jayenne; git pull"
else
  run "cd ${REGDIR}; git clone git@gitlab.lanl.gov:jayenne/jayenne.git"
fi
# Capsaicin
echo " "
echo "Updating $REGDIR/capsaicin..."
if test -d ${REGDIR}/capsaicin/scripts; then
run "cd ${REGDIR}/capsaicin/scripts; ${SVN} update"
else
  run "mkdir -p ${REGDIR}/capsaicin; cd ${REGDIR}/capsaicin"
  run "${SVN} co svn+ssh://ccscs7.lanl.gov/ccs/codes/radtran/svn/capsaicin/trunk/scripts"
fi

##---------------------------------------------------------------------------##
## End update_regression_scripts.sh
##---------------------------------------------------------------------------##
