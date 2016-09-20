#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/sync_repository.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# This script is used for 2 similar but distinct operations:
#
# 1. It mirrors portions of the Capsaicin SVN repositories to a HPC
#    location. The repository must be mirrored because the ctest regression
#    system must be run from the HPC backend (via msub) where access to
#    ccscs7:/ccs/codes/radtran/svn is not available.
# 2. It also mirrors git@github.com/losalamos/Draco.git and
#    git@gitlab.lanl.gov/jayenne/jayenne.git to these locations:
#    - ccscs7:/ccs/codes/radtran/git
#    - darwin-fe:/usr/projects/draco/regress/git
#    On ccscs7, this is done to allow Redmine to parse the current repository
#    preseting a GUI view and scraping commit information that connects to
#    tracked issues. On darwin, this is done to allow the regressions running on
#    the compute node to access the latest git repository. This also copies down
#    all pull requests.

target="`uname -n | sed -e s/[.].*//`"

# Locate the directory that this script is located in:
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#
# MODULES
#
# Determine if the module command is available
modcmd=`declare -f module`
# If not found, look for it in /usr/share/Modules (ML)
if [[ ! ${modcmd} ]]; then
  if test -f /usr/share/Modules/init/bash; then
    source /usr/share/Modules/init/bash
  else
    echo "ERROR: The module command was not found. No modules will be loaded."
  fi
fi
modcmd=`declare -f module`

#
# Environment
#

# import some bash functions
source $scriptdir/scripts/common.sh

# Ensure that the permissions are correct
run "umask 0002"
svnhostmachine=ccscs7

case ${target} in
  ccscs*)
    run "module load user_contrib subversion git"
    regdir=/scratch/regress
    gitroot=/ccs/codes/radtran/git
    VENDOR_DIR=/scratch/vendors
    keychain=keychain-2.8.2
    ;;
  darwin-fe* | cn[0-9]*)
    regdir=/usr/projects/draco/regress
    gitroot=$regdir/git
    VENDOR_DIR=/usr/projects/draco/vendors
    keychain=keychain-2.7.1
    ;;
  *)
    # HPC - Moonlight.
    run "module load user_contrib svn git"
    regdir=/usr/projects/jayenne/regress
    gitroot=$regdir/git
    svnroot=$regdir/svn
    VENDOR_DIR=/usr/projects/draco/vendors
    keychain=keychain-2.7.1
    ;;
esac

if ! test -d $regdir; then
  mkdir -p $regdir
fi

# Credentials via Keychain (SSH)
# http://www.cyberciti.biz/faq/ssh-passwordless-login-with-keychain-for-scripts
MYHOSTNAME="`uname -n`"
$VENDOR_DIR/$keychain/keychain $HOME/.ssh/cmake_dsa
if test -f $HOME/.keychain/$MYHOSTNAME-sh; then
  run "source $HOME/.keychain/$MYHOSTNAME-sh"
else
  echo "Error: could not find $HOME/.keychain/$MYHOSTNAME-sh"
fi

# ---------------------------------------------------------------------------- #
# Create copies of SVN and GIT repositories on the local file system
# ---------------------------------------------------------------------------- #

case ${target} in
  ccscs*)
    # Keep local (ccscs7:/ccs/codes/radtran/git) copies of the github and gitlab
    # repositories. This location can be parsed by redmine.

    echo " "
    echo "Copy Draco git repository to the local file system..."
    if test -d $gitroot/Draco.git; then
      run "cd $gitroot/Draco.git"
      run "git fetch origin +refs/heads/*:refs/heads/*"
      run "git fetch origin +refs/pull/*:refs/pull/*"
      run "git reset --soft"
    else
      run "mkdir -p $gitroot"
      run "cd $gitroot"
      run "git clone --bare git@github.com:losalamos/Draco.git Draco.git"
    fi

    echo " "
    echo "Copy Jayenne git repository to the local file system..."
    if test -d $gitroot/jayenne.git; then
      run "cd $gitroot/jayenne.git"
      run "git fetch origin +refs/heads/*:refs/heads/*"
      run "git fetch origin +refs/merge-requests/*:refs/merge-requests/*"
      run "git reset --soft"
    else
      run "mkdir -p $gitroot; cd $gitroot"
      run "git clone --bare git@gitlab.lanl.gov:jayenne/jayenne.git jayenne.git"
    fi
    ;;
  darwin-fe* | ml-fey*)
    # Keep local (ccscs7:/ccs/codes/radtran/git) copies of the github and gitlab
    # repositories. This location can be parsed by redmine. For darwin, the
    # backend can't see gitlab, so keep a copy of the repository local.

    echo " "
    echo "Copy Draco git repository to the local file system..."
    if test -d $gitroot/Draco.git; then
      run "cd $gitroot/Draco.git"
      run "git fetch origin +refs/heads/*:refs/heads/*"
      run "git fetch origin +refs/pull/*:refs/pull/*"
      run "git reset --soft"
    else
      run "mkdir -p $gitroot"
      run "cd $gitroot"
      run "git clone --bare git@github.com:losalamos/Draco.git Draco.git"
    fi

    echo " "
    echo "Copy Jayenne git repository to the local file system..."
    if test -d $gitroot/jayenne/jayenne.git; then
      run "cd $gitroot/jayenne/jayenne.git"
      run "git fetch origin +refs/heads/*:refs/heads/*"
      run "git fetch origin +refs/merge-requests/*:refs/merge-requests/*"
      run "git reset --soft"
    else
      run "mkdir -p $gitroot/jayenne; cd $gitroot/jayenne"
      run "git clone --bare git@gitlab.lanl.gov:jayenne/jayenne.git jayenne.git"
    fi
    ;;
  *)
    #
    # HPC: Mirror the capsaicin svn repository
    #
    if ! test -d $svnroot; then
      echo "*** ERROR ***"
      echo "*** SVN repository not found ***"
      exit 1
      # http://journal.paul.querna.org/articles/2006/09/14/using-svnsync/
      # mkdir -p ${svnroot}; cd ${svnroot}
      # svnadmin create ${svnroot}/jayenne
      # chgrp -R draco jayenne; chmod -R g+rwX,o=g-w jayenne
      # cd jayenne/hooks
      # cp pre-commit.tmpl pre-commit; chmod 775 pre-commit
      # vi pre-commit; comment out all code and add...
      #if ! test `whoami` = 'kellyt'; then
      #echo "This is a read only repository.  The real SVN repository is"
      #echo "at svn+ssh://ccscs8/ccs/codes/radtran/svn/draco."
      #exit 1
      #fi
      #exit 0
      # cp pre-revprop-change.tmpl pre-revprop-change; chmod 775 \
      #    pre-revprop-change
      # vi pre-revprop-change --> comment out all code.
      # cd $svnroot
      # svnsync init file:///${svnroot}/jayenne svn+ssh://ccscs8/ccs/codes/radtran/svn/jayenne
      # svnsync sync file:///${svnroot}/jayenne
    fi

    echo " "
    echo "Copy capsaicin svn repository to the local file system..."
    run "svnsync --non-interactive sync file:///${svnroot}/capsaicin"
    ;;
esac

#------------------------------------------------------------------------------#
# End sync_repository.sh
#------------------------------------------------------------------------------#
