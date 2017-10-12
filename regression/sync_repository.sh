#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/sync_repository.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# switch to group 'ccsrad' and set umask
if [[ $(id -gn) != ccsrad ]]; then
  exec sg ccsrad "$0 $*"
fi
umask 0007

# Locate the directory that this script is located in:
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
target="`uname -n | sed -e s/[.].*//`"

# Prevent multiple copies of this script from running at the same time:
lockfile=/var/tmp/sync_repository_$target.lock
[ "${FLOCKER}" != "${lockfile}" ] && exec env FLOCKER="${lockfile}" flock -en "${lockfile}" "${0}" "$@" || :


# All output will be saved to this log file.  This is also the lockfile for flock.
timestamp=`date +%Y%m%d-%H%M`
logdir="$( cd $scriptdir/../../logs && pwd )"
logfile=$logdir/sync_repository_${target}_${timestamp}.log

# Debug stuff
verbose=off
if [[ ${verbose:-off} == "on" ]]; then
  echo "looking for locks..."
  echo "   FLOCKER     = ${FLOCKER}"
  echo "   lockfile    = $lockfile"
  echo "   logfile     = $logfile"
  echo "   script name = ${0}"
  echo "   script args = $@"
fi

# Redirect all future output to the logfile.
exec > $logfile
exec 2>&1

# import some bash functions
source $scriptdir/scripts/common.sh

echo -e "Executing $0 $*...\n"
echo "Group: `id -gn`"
echo -e "umask: `umask`\n"

#------------------------------------------------------------------------------#
# This script is used for 2 similar but distinct operations:
#
# 1. It mirrors git@github.com/lanl/Draco.git,
#    git@gitlab.lanl.gov/jayenne/jayenne.git and
#    git@gitlab.lanl.gov/capsaicin/capsaicin to these locations:
#    - ccscs7:/ccs/codes/radtran/git
#    - darwin-fe:/usr/projects/draco/regress/git
#    On ccscs7, this is done to allow Redmine to parse the current repository
#    preseting a GUI view and scraping commit information that connects to
#    tracked issues. On darwin, this is done to allow the regressions running on
#    the compute node to access the latest git repository. This also copies down
#    all pull requests.
# 2. It captures the output produced during the mirroring process.  If a new PR
#    is found in the mirrored repository, continuous integration testing is
#    started.

#
# MODULES
#

# If not found, look for it in /usr/share/Modules (ML)
if [[ `fn_exists module` == 0 ]]; then
  case ${target} in
    tt-fey*) module_init_dir=/opt/cray/pe/modules/3.2.10.4/init/bash ;;
    # snow (Toss3, lmod)
    sn-fey*) module_init_dir=/usr/share/lmod/lmod/init/profile ;;
    # ccs-net (lmod)
    ccscs*)  module_init_dir=/usr/share/lmod/lmod/init/bash ;;
    # darwin, ml
    *)       module_init_dir=/usr/share/Modules/init/bash ;;
  esac
  if [[ -f ${module_init_dir} ]]; then
    echo "Module environment not found, trying to init the module environment."
    run "source ${module_init_dir}"
  else
    echo "ERROR: The module command was not found. No modules will be loaded."
  fi
  if [[ `fn_exists module` == 0 ]]; then
    echo "ERROR: the module command was not found (even after sourcing $module_init_dir"
    exit 1
  fi
fi

#
# Environment
#

case ${target} in
  ccscs*)
    run "module use /usr/share/lmod/lmod/modulefiles/Core"
    run "module load user_contrib subversion git"
    regdir=/scratch/regress
    gitroot=/ccs/codes/radtran/git.${target}
    VENDOR_DIR=/scratch/vendors
    keychain=keychain-2.8.2
    ;;
  darwin-fe* | cn[0-9]*)
    regdir=/usr/projects/draco/regress
    gitroot=$regdir/git
    VENDOR_DIR=/usr/projects/draco/vendors
    keychain=keychain-2.7.1
    ;;
  ml-fey*)
    run "module load user_contrib subversion git"
    regdir=/usr/projects/jayenne/regress
    gitroot=$regdir/git
    VENDOR_DIR=/usr/projects/draco/vendors
    keychain=keychain-2.7.1
    ;;
  sn-fey*)
    run "module use --append /usr/projects/hpcsoft/toss3/modulefiles/snow/compiler"
    run "module use --append /usr/projects/hpcsoft/toss3/modulefiles/snow/libraries"
    run "module use --append /usr/projects/hpcsoft/toss3/modulefiles/snow/misc"
    run "module use --append /usr/projects/hpcsoft/toss3/modulefiles/snow/mpi"
    run "module use --append /usr/projects/hpcsoft/toss3/modulefiles/snow/tools"
    run "module load user_contrib subversion git"
    regdir=/usr/projects/jayenne/regress
    gitroot=$regdir/git.sn
    VENDOR_DIR=/usr/projects/draco/vendors
    keychain=keychain-2.7.1
    ;;
  tt-fey*)
    run "module use /usr/projects/hpcsoft/cle6.0/modulefiles/trinitite/misc"
    run "module use /usr/projects/hpcsoft/cle6.0/modulefiles/trinitite/tools"
    run "module load user_contrib subversion git"
    regdir=/usr/projects/jayenne/regress
    gitroot=$regdir/git.tt
    VENDOR_DIR=/usr/projects/draco/vendors
    keychain=keychain-2.7.1
    ;;
esac

if ! [[ -d $regdir ]]; then
  run "mkdir -p $regdir"
fi

# Credentials via Keychain (SSH)
# http://www.cyberciti.biz/faq/ssh-passwordless-login-with-keychain-for-scripts
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
# Create copies of GIT repositories on the local file system
#
# Locations:
# - ccs-net servers:
#   /ccs/codes/radtran/git
# - HPC (moonlight, snow, trinitite, etc.)
#   /usr/projects/draco/jayenne/regress/git
#
# Keep local copies of the github, gitlab and svn repositories. This local
# filesystem location is needed for our regression system on some systems where
# the back-end worker nodes can't see the outside world. Additionally, the
# output produced by downloading the repository to the local filesystem can be
# parsed so that CI regressions can be started. On the CCS-NET, the local copy
# is also visible to Redmine so changes in the repository will show up in the
# datase (GUI repository, wiki/ticket references to commits).
# ---------------------------------------------------------------------------- #

# DRACO: For all machines running this script, copy all of the git repositories
# to the local file system.

# Store some output into a local file to simplify parsing.
TMPFILE_DRACO_HEAD=$(mktemp /var/tmp/draco_repo_sync.XXXXXXXXXX) || die "Failed to create temporary file"
TMPFILE_DRACO_PULL=$(mktemp /var/tmp/draco_repo_sync.XXXXXXXXXX) || die "Failed to create temporary file"

echo " "
echo "Copy Draco git repository to the local file system..."
if [[ -d $gitroot/Draco.git ]]; then
  run "cd $gitroot/Draco.git"
  run "git fetch origin +refs/heads/*:refs/heads/*" &> $TMPFILE_DRACO_HEAD
  run "cat $TMPFILE_DRACO_HEAD"
  run "git fetch origin +refs/pull/*:refs/pull/*" &> $TMPFILE_DRACO_PULL
  run "cat $TMPFILE_DRACO_PULL"
  run "git reset --soft"
  run "chgrp -R draco $gitroot/Draco.git"
  run "chmod -R g+rwX $gitroot/Draco.git"
else
  run "mkdir -p $gitroot"
  run "cd $gitroot"
  run "git clone --bare git@github.com:lanl/Draco.git Draco.git"
  run "cd Draco.git"
  run "git fetch origin +refs/heads/*:refs/heads/*"
  run "git fetch origin +refs/pull/*:refs/pull/*"
fi

# JAYENNE: For all machines running this scirpt, copy all of the git repositories
# to the local file system.

# Store some output into a local file to simplify parsing.
TMPFILE_JAYENNE=$(mktemp /var/tmp/jayenne_repo_sync.XXXXXXXXXX) || die "Failed to create temporary file"
echo " "
echo "Copy Jayenne git repository to the local file system..."
if [[ -d $gitroot/jayenne.git ]]; then
  run "cd $gitroot/jayenne.git"
  run "git fetch origin +refs/heads/*:refs/heads/*"
  run "git fetch origin +refs/merge-requests/*:refs/merge-requests/*" &> $TMPFILE_JAYENNE
  run "cat $TMPFILE_JAYENNE"
  run "git reset --soft"
  run "chgrp ccsrad $gitroot/jayenne.git"
  run "chmod -R g+rwX,o-rwX $gitroot/jayenne.git"
else
  run "mkdir -p $gitroot; cd $gitroot"
  run "git clone --bare git@gitlab.lanl.gov:jayenne/jayenne.git jayenne.git"
  run "cd jayenne.git"
  run "git fetch origin +refs/heads/*:refs/heads/*"
  run "git fetch origin +refs/merge-requests/*:refs/merge-requests/*"
fi

# CAPSAICIN: For all machines running this scirpt, copy all of the git repositories
# to the local file system.

# Store some output into a local file to simplify parsing.
TMPFILE_CAPSAICIN=$(mktemp /var/tmp/capsaicin_repo_sync.XXXXXXXXXX) || die "Failed to create temporary file"
echo " "
echo "Copy Capsaicin git repository to the local file system..."
if [[ -d $gitroot/capsaicin.git ]]; then
  run "cd $gitroot/capsaicin.git"
  run "git fetch origin +refs/heads/*:refs/heads/*"
  run "git fetch origin +refs/merge-requests/*:refs/merge-requests/*" &> $TMPFILE_CAPSAICIN
  run "cat $TMPFILE_CAPSAICIN"
  run "git reset --soft"
  run "chgrp ccsrad $gitroot/capsaicin.git"
else
  run "mkdir -p $gitroot; cd $gitroot"
  run "git clone --bare git@gitlab.lanl.gov:capsaicin/capsaicin.git capsaicin.git"
  run "cd capsaicin.git"
  run "git fetch origin +refs/heads/*:refs/heads/*"
  run "git fetch origin +refs/merge-requests/*:refs/merge-requests/*"
fi

#------------------------------------------------------------------------------#
# Mirror git repository for redmine integration
#------------------------------------------------------------------------------#

if [[ ${target} == "ccscs7" ]]; then

  # Keep a copy of the bare repo for Redmine.  This version doesn't have the
  # PRs since this seems to confuse Redmine.
  projects="jayenne capsaicin"
  for p in $projects; do
    echo -e "\n(Redmine) Copy ${p} git repository to the local file system..."
    if [[ -d $gitroot/${p}-redmine.git ]]; then
      run "cd $gitroot/${p}-redmine.git"
      run "git fetch origin +refs/heads/*:refs/heads/*"
      run "git reset --soft"
      run "cd $gitroot"
      run "chmod -R g+rwX,o-rwX $gitroot/${p}-redmine.git"
    else
      run "mkdir -p $gitroot"
      run "cd $gitroot"
      if [[ ${p} == "Draco" ]]; then
        run "git clone --bare git@github.com:lanl/${p}.git ${p}-redmine.git"
      else
        run "git clone --bare git@gitlab.lanl.gov:${p}/${p}.git ${p}-redmine.git"
      fi
      run "chmod -R g+rwX,o-rwX ${p}-redmine.git"
      run "chmod g+s ${p}-redmine.git"
    fi
  done

fi

#------------------------------------------------------------------------------#
# Continuous Integration Hooks:
#
# Extract a list of PRs that are new and optionally start regression run by
# parsing the output of 'git fetch' from above.
#
# The output from the above command may include text of the form:
#
#   [new ref]        refs/pull/84/head -> refs/pull/84/head <-- new PR
#   03392b8..fd3eabc refs/pull/86/head -> refs/pull/86/head <-- updated PR
#   881a1f4...86c80c8 refs/pull/157/merge -> refs/pull/157/merge  (forced update)
#
# From gitlab:
#
# * [new ref]        refs/merge-requests/141/head -> refs/merge-requests/141/head
#
# Extract a list of PRs that are new and optionally start regression run
# ------------------------------------------------------------------------------#

echo " "
echo "========================================================================"
echo "Starting CI regressions (if any)"
date
echo "========================================================================"
echo " "
# Draco CI ------------------------------------------------------------

# Did we find 'merge' in the repo sync?  If so, then we should reset the
# last-draco tagfile.
unset rmlastdracotag
draco_merge=`cat $TMPFILE_DRACO_HEAD | grep -c "develop    -> develop"`
if [[ $draco_merge -gt 0 ]]; then
  rmlastdracotag="-t"
fi
draco_prs=`cat $TMPFILE_DRACO_PULL | grep -e 'refs/pull/[0-9]*/\(head\|merge\)' | sed -e 's%.*/\([0-9][0-9]*\)/.*%\1%'`
# remove any duplicates
draco_prs=`echo $draco_prs | xargs -n1 | sort -u | xargs`
for pr in $draco_prs; do
  run "$scriptdir/checkpr.sh -r -p draco -f $pr $rmlastdracotag"
done

# Jayenne CI ----------------------------------------------------------

jayenne_prs=`cat $TMPFILE_JAYENNE | grep -e 'refs/merge-requests/[0-9]*/\(head\|merge\)' | sed -e 's%.*/\([0-9][0-9]*\)/.*%\1%'`
# remove any duplicates
jayenne_prs=`echo $jayenne_prs | xargs -n1 | sort -u | xargs`
for pr in $jayenne_prs; do
  run "$scriptdir/checkpr.sh -r -p jayenne -f $pr"
done

# Capsaicin CI ----------------------------------------------------------

capsaicin_prs=`cat $TMPFILE_CAPSAICIN | grep -e 'refs/merge-requests/[0-9]*/\(head\|merge\)' | sed -e 's%.*/\([0-9][0-9]*\)/.*%\1%'`
# remove any duplicates
capsaicin_prs=`echo $capsaicin_prs | xargs -n1 | sort -u | xargs`
for pr in $capsaicin_prs; do
  run "$scriptdir/checkpr.sh -r -p capsaicin -f $pr"
done

# Wait for all subprocesses to finish before exiting this script
if [[ `jobs -p | wc -l` -gt 0 ]]; then
  echo " "
  echo "Jobs still running (if any):"
  for job in `jobs -p`; do
    echo "  waiting for job $job to finish..."
    wait $job
    echo "  waiting for job $job to finish...done"
  done
fi

# Cleanup
echo " "
echo "Cleaning up..."
run "rm $TMPFILE_DRACO_PULL $TMPFILE_DRACO_HEAD $TMPFILE_JAYENNE $TMPFILE_CAPSAICIN"
run "rm $lockfile"

echo " "
date
echo "All done."

#------------------------------------------------------------------------------#
# End sync_repository.sh
#------------------------------------------------------------------------------#
