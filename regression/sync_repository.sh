#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/sync_repository.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2019, Triad National Security, LLC.
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
if ! [[ -d /var/tmp/$USER ]]; then
  mkdir -p /var/tmp/$USER || exit 2
fi
lockfile=/var/tmp/$USER/sync_repository_$target.lock
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
#    git@gitlab.lanl.gov/capsaicin/core to these locations:
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
    # darwin
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
    run "module load user_contrib git"
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
  sn-fey*)
    run "module use --append /usr/projects/hpcsoft/toss3/modulefiles/snow/compiler"
    run "module use --append /usr/projects/hpcsoft/toss3/modulefiles/snow/libraries"
    run "module use --append /usr/projects/hpcsoft/toss3/modulefiles/snow/misc"
    run "module use --append /usr/projects/hpcsoft/toss3/modulefiles/snow/mpi"
    run "module use --append /usr/projects/hpcsoft/toss3/modulefiles/snow/tools"
    run "module load user_contrib git"
    regdir=/usr/projects/jayenne/regress
    gitroot=$regdir/git.sn
    VENDOR_DIR=/usr/projects/draco/vendors
    keychain=keychain-2.7.1
    ;;
  tt-fey*)
    run "module use /usr/projects/hpcsoft/cle6.0/modulefiles/trinitite/misc"
    run "module use /usr/projects/hpcsoft/cle6.0/modulefiles/trinitite/tools"
    run "module load user_contrib git"
    regdir=/usr/projects/jayenne/regress
    gitroot=$regdir/git.tt
    VENDOR_DIR=/usr/projects/draco/vendors
    keychain=keychain-2.8.2
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
# - HPC (snow, trinitite, etc.)
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

# associative array to keep track of the TMPFILES
declare -A tmpfiles=()

# List of repositories (also used by push_repositories_xf.sh and
# pull_repositories_xf.sh).  It defines $github_projects and $gitlab_projects.
source ${scriptdir}/repository_list.sh

# Github.com/lanl repositories:

for project in ${github_projects[@]}; do

  namespace=`echo $project | sed -e 's%/.*%%'`
  repo=`echo $project | sed -e 's%.*/%%'`

  # Store some output into a local file to simplify parsing.
  tmpfiles[${project}]=$(mktemp /var/tmp/$USER/${namespace}_${repo}_repo_sync.XXXXXXXXXX) || die "Failed to create temporary file"

  echo -e "\nCopy ${project}'s git repository to the local file system...\n"
  if [[ -d $gitroot/${project}.git ]]; then
    run "cd $gitroot/${project}.git"
    run "git fetch origin +refs/heads/*:refs/heads/* &> ${tmpfiles[${project}]}"
    run "git fetch origin +refs/pull/*:refs/pull/* >> ${tmpfiles[${project}]} 2>&1"
    run "cat ${tmpfiles[${project}]}"
    run "git reset --soft"
  else
    run "mkdir -p $gitroot/$namespace; cd $gitroot/$namespace"
    run "git clone --mirror git@github.com:${project}.git ${repo}.git"
    run "cd ${repo}.git"
    run "git fetch origin +refs/heads/*:refs/heads/*"
    run "git fetch origin +refs/pull/*:refs/pull/*"
    # if this is a brand new checkout, then do not run any ci (there might be hundreds).
    no_ci=yes
  fi
  run "chgrp -R draco $gitroot/${namespace}"
  run "chmod -R g+rwX,o=g-w $gitroot/${namespace}"
  run "find draco -type d -exec chmod g+s {} \;"
done

# Gitlab.lanl.gov repositories:

for project in ${gitlab_projects[@]}; do

  namespace=`echo $project | sed -e 's%/.*%%'`
  repo=`echo $project | sed -e 's%.*/%%'`

  # Store some output into a local file to simplify parsing.
  tmpfiles[${project}]=$(mktemp /var/tmp/$USER/${namespace}_${repo}_repo_sync.XXXXXXXXXX) || die "Failed to create temporary file"

  echo -e "\nCopy ${project}'s git repository to the local file system...\n"
  if [[ -d $gitroot/${project}.git ]]; then
    run "cd $gitroot/${project}.git"
    run "git fetch origin +refs/heads/*:refs/heads/*"
    run "git fetch origin +refs/merge-requests/*:refs/merge-requests/* &> ${tmpfiles[${project}]}"
    run "cat ${tmpfiles[${project}]}"
    run "git reset --soft"
  else
    run "mkdir -p $gitroot/$namespace; cd $gitroot/$namespace"
    run "git clone --mirror git@gitlab.lanl.gov:${project}.git ${repo}.git"
    run "cd ${repo}.git"
    run "git fetch origin +refs/heads/*:refs/heads/*"
    run "git fetch origin +refs/merge-requests/*:refs/merge-requests/*"
    # if this is a brand new checkout, then do not run any ci (there might be hundreds).
    no_ci=yes
  fi
  run "chgrp -R ccsrad $gitroot/${namespace}"
  run "chmod -R g+rwX,o-rwX $gitroot/${namespace}"
  run "find $gitroot/$namespace -type d -exec chmod g+s {} \;"

done

#------------------------------------------------------------------------------#
# Mirror git repository for redmine integration
#------------------------------------------------------------------------------#

if [[ ${target} == "ccscs7" ]]; then

  # Keep a copy of the bare repo for Redmine.  This version doesn't have the
  # PRs since this seems to confuse Redmine.
  redmine_projects="jayenne core trt npt"
  for p in $redmine_projects; do
    echo -e "\n(Redmine) Copy ${p} git repository to the local file system..."
    if [[ -d $gitroot/${p}-redmine.git ]]; then
      run "cd $gitroot/${p}-redmine.git"
      run "git fetch origin +refs/heads/*:refs/heads/*"
      run "git reset --soft"
      run "cd $gitroot"
      run "chmod -R g+rwX,o-rwX $gitroot/${p}-redmine.git"
      run "find $gitroot/${p}-redmine.git -type d -exec chmod g+s {} \;"

    else
      run "mkdir -p $gitroot"
      run "cd $gitroot"
      if [[ ${p} == "Draco" ]]; then
        run "git clone --bare git@github.com:lanl/${p}.git ${p}-redmine.git"
      else
        run "git clone --bare git@gitlab.lanl.gov:${p}/${p}.git ${p}-redmine.git"
      fi
      run "chmod -R g+rwX,o-rwX ${p}-redmine.git"
      run "find ${p}-redmine.git -type d -exec chmod g+s {} \;"
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

echo -e "\n========================================================================"
echo -e "Starting CI regressions (if any)"
date
echo -e "========================================================================"

# if this is a brand new checkout, then do not run any ci (there might be hundreds).
if [[ ${no_ci:-no} == yes ]]; then
  dry_run=yes
  echo "Enable dry_run=yes mode."
fi

# CI ------------------------------------------------------------

for project in ${git_projects[@]}; do

  namespace=`echo $project | sed -e 's%/.*%%'`
  repo=`echo $project | sed -e 's%.*/%%'`

  # Did we find 'merge' in the repo sync?  If so, then we should reset the
  # last-draco tagfile.
  unset rm_last_build

  echo -e "\n----- ${project} -----\n"

  case $repo in
    Draco)
      if [[ `cat ${tmpfiles[${project}]} | grep -c "develop    -> develop"` -gt 0 ]]; then
        rm_last_build="-t"
      fi
      # Extract PR number (if any)
      prs=`cat ${tmpfiles[${project}]} | grep -e 'refs/pull/[0-9]*/\(head\|merge\)' | sed -e 's%.*/\([0-9][0-9]*\)/.*%\1%'`
      ;;

    jayenne | core | trt | npt)
      # Extract PR number (if any)
      prs=`cat ${tmpfiles[${project}]} | grep -e 'refs/merge-requests/[0-9]*/\(head\|merge\)' | sed -e 's%.*/\([0-9][0-9]*\)/.*%\1%'`
      ;;
  esac

  case $repo in
    Draco | jayenne | core | trt | npt)
      repo_lc=`echo $repo | tr '[:upper:]' '[:lower:]'`
      # remove any duplicates
      prs=`echo $prs | xargs -n1 | sort -u | xargs`
      for pr in $prs; do
        run "$scriptdir/checkpr.sh -r -p ${repo_lc} -f $pr $rm_last_build"
      done
      ;;
    *)
      echo " - No CI defined."
      ;;
  esac

done

#------------------------------------------------------------------------------#
# Wait, cleanup, done
#------------------------------------------------------------------------------#

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

unset dry_run

for file in "${tmpfiles[@]}"; do
  run "rm $file"
done
run "rm $lockfile"

echo " "
date
echo "All done."

#------------------------------------------------------------------------------#
# End sync_repository.sh
#------------------------------------------------------------------------------#
