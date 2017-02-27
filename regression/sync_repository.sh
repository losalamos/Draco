#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/sync_repository.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# This script is used for 2 similar but distinct operations:
#
# 1. It mirrors portions of the Capsaicin SVN repositories to a HPC
#    location. The repository must be mirrored because the ctest regression
#    system must be run from the HPC backend (via msub) where access to
#    ccscs7:/ccs/codes/radtran/svn is not available.
# 2. It also mirrors git@github.com/lanl/Draco.git and
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

# All output will be saved to this log file.  This is also the lockfile for flock.
logdir="$( cd $scriptdir/../../logs && pwd )"
logfile=$logdir/sync_repository_$target.log
lockfile=/var/tmp/sync_repository_$target.lock

# Prevent multiple copies of this script from running at the same time:
[ "${FLOCKER}" != "${lockfile}" ] && exec env FLOCKER="${lockfile}" flock -en "${lockfile}" "${0}" "$@" || :

# Redirect all future output to the logfile.
exec > $logfile
exec 2>&1

#
# MODULES
#
# Determine if the module command is available
modcmd=`declare -f module`
# If not found, look for it in /usr/share/Modules (ML)
if [[ ! ${modcmd} ]]; then
  case ${target} in
    tt-fey*) module_init_dir=/opt/cray/pe/modules/3.2.10.4/init/bash ;;
    # snow (Toss3)
    sn-fey*) module_init_dir=/usr/share/lmod/lmod/init/profile ;;
    # ccs-net, darwin, ml
    *)       module_init_dir=/usr/share/Modules/init/bash ;;
  esac
  if test -f ${module_init_dir}; then
    source ${module_init_dir}
  else
    echo "ERROR: The module command was not found. No modules will be loaded."
  fi
  modcmd=`declare -f module`
  if [[ ! ${modcmd} ]]; then
    echo "ERROR: the module command was not found (even after sourcing $module_init_dir"
    exit 1
  fi
fi

#
# Environment
#

##----------------------------------------------------------------------------##
# Pause until the 'last modified' timestamp of file $1 to be $2 seconds old.
function allow_file_to_age
{
  if [[ ! $2 ]]; then
    echo "ERROR: This function requires two arguments: a filename and an age value (sec)."
    exit 1
  fi

  # If file does not exist, no need to wait.
  if [[ ! -f $1 ]]; then
    return
  fi

  # assume file was last modified 0 seconds ago.
  timediff=0

  # If no changes for $2 seconds, continue
  # else, wait until until file, $1, hasn't been touched for $2 seconds.
  while [[ $timediff -lt $2 ]]; do
    eval "$(date +'now=%s')"
    pr_last_check=$(date +%s -r $1)
    timediff=$(expr $now - $pr_last_check)
    sleep 30s
  done
}

# import some bash functions
source $scriptdir/scripts/common.sh

# Ensure that the permissions are correct
run "umask 0002"

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

# DRACO: For all machines running this scirpt, copy all of the git repositories
# to the local file system.

# Store some output into a local file to simplify parsing.
TMPFILE_DRACO=$(mktemp /var/tmp/draco_repo_sync.XXXXXXXXXX) || { echo "Failed to create temporary file"; exit 1; }

echo " "
echo "Copy Draco git repository to the local file system..."
if test -d $gitroot/Draco.git; then
  run "cd $gitroot/Draco.git"
  run "git fetch origin +refs/heads/*:refs/heads/*"
  run "git fetch origin +refs/pull/*:refs/pull/*" &> $TMPFILE_DRACO
  cat $TMPFILE_DRACO
  run "git reset --soft"
else
  run "mkdir -p $gitroot"
  run "cd $gitroot"
  run "git clone --bare git@github.com:lanl/Draco.git Draco.git"
  run "Draco.git"
  run "git fetch origin +refs/heads/*:refs/heads/*"
  run "git fetch origin +refs/pull/*:refs/pull/*"
fi
case ${target} in
  ccscs7*)
    # Keep a copy of the bare repo for Redmine.  This version doesn't have the
    # PRs since this seems to confuse Redmine.
    echo " "
    echo "(Redmine) Copy Draco git repository to the local file system..."
    if test -d $gitroot/Draco-redmine.git; then
      run "cd $gitroot/Draco-redmine.git"
      run "git fetch origin +refs/heads/*:refs/heads/*"
      run "git reset --soft"
    else
      run "mkdir -p $gitroot"
      run "cd $gitroot"
      run "git clone --mirror git@github.com:lanl/Draco.git Draco-redmine.git"
      run "chmod -R g+rwX Draco-redmine.git"
    fi
    ;;
esac


# JAYENNE: For all machines running this scirpt, copy all of the git repositories
# to the local file system.

# Store some output into a local file to simplify parsing.
TMPFILE_JAYENNE=$(mktemp /var/tmp/jayenne_repo_sync.XXXXXXXXXX) || { echo "Failed to create temporary file"; exit 1; }
echo " "
echo "Copy Jayenne git repository to the local file system..."
if test -d $gitroot/jayenne.git; then
  run "cd $gitroot/jayenne.git"
  run "git fetch origin +refs/heads/*:refs/heads/*"
  run "git fetch origin +refs/merge-requests/*:refs/merge-requests/*" &> $TMPFILE_JAYENNE
  cat $TMPFILE_JAYENNE
  run "git reset --soft"
else
  run "mkdir -p $gitroot; cd $gitroot"
  run "git clone --bare git@gitlab.lanl.gov:jayenne/jayenne.git jayenne.git"
  run "cd jayenne.git"
  run "git fetch origin +refs/heads/*:refs/heads/*"
  run "git fetch origin +refs/merge-requests/*:refs/merge-requests/*"
fi
case ${target} in
  ccscs7*)
    # Keep a copy of the bare repo for Redmine.  This version doesn't have the
    # PRs since this seems to confuse Redmine.
    echo " "
    echo "(Redmine) Copy Jayenne git repository to the local file system..."
    if test -d $gitroot/jayenne-redmine.git; then
      run "cd $gitroot/jayenne-redmine.git"
      run "git fetch origin +refs/heads/*:refs/heads/*"
      run "git reset --soft"
    else
      run "mkdir -p $gitroot"
      run "cd $gitroot"
      run "git clone --mirror git@gitlab.lanl.gov:jayenne/jayenne.git jayenne-redmine.git"
      run "chmod -R g+rwX jayenne-redmine.git"
    fi
    ;;
esac

# CAPSAICIN: For all machines running this scirpt, copy all of the git repositories
# to the local file system.

# Store some output into a local file to simplify parsing.
TMPFILE_CAPSAICIN=$(mktemp /var/tmp/capsaicin_repo_sync.XXXXXXXXXX) || { echo "Failed to create temporary file"; exit 1; }
echo " "
echo "Copy Capsaicin git repository to the local file system..."
if test -d $gitroot/capsaicin.git; then
  run "cd $gitroot/capsaicin.git"
  run "git fetch origin +refs/heads/*:refs/heads/*"
  run "git fetch origin +refs/merge-requests/*:refs/merge-requests/*" &> $TMPFILE_CAPSAICIN
  cat $TMPFILE_CAPSAICIN
  run "git reset --soft"
else
  run "mkdir -p $gitroot; cd $gitroot"
  run "git clone --bare git@gitlab.lanl.gov:capsaicin/capsaicin.git capsaicin.git"
  run "cd capsaicin.git"
  run "git fetch origin +refs/heads/*:refs/heads/*"
  run "git fetch origin +refs/merge-requests/*:refs/merge-requests/*"
fi
case ${target} in
  ccscs7*)
    # Keep a copy of the bare repo for Redmine.  This version doesn't have the
    # PRs since this seems to confuse Redmine.
    echo " "
    echo "(Redmine) Copy Capsaicin git repository to the local file system..."
    if test -d $gitroot/capsaicin-redmine.git; then
      run "cd $gitroot/capsaicin-redmine.git"
      run "git fetch origin +refs/heads/*:refs/heads/*"
      run "git reset --soft"
    else
      run "mkdir -p $gitroot"
      run "cd $gitroot"
      run "git clone --mirror git@gitlab.lanl.gov:capsaicin/capsaicin.git capsaicin-redmine.git"
      run "chmod -R g+rwX capsaicin-redmine.git"
    fi
    ;;
esac

#------------------------------------------------------------------------------#
# Continuous Integration Hooks:
#
# Extract a list of PRs that are new and optionally start regression run by
# parsing the output of 'git fetch' from above.
#
# The output from the above command may include text of the form:
#   [new ref]        refs/pull/84/head -> refs/pull/84/head <-- new PR
#   03392b8..fd3eabc refs/pull/86/head -> refs/pull/86/head <-- updated PR
#                    Extract a list of PRs that are new and optionally start
#                    regression run
# ------------------------------------------------------------------------------#

# Draco CI ------------------------------------------------------------
draco_prs=`grep 'refs/pull/[0-9]*/head$' $TMPFILE_DRACO | awk '{print $NF}'`
for prline in $draco_prs; do
  echo " "
  case ${target} in

    # CCS-NET: Coverage (Debug) & Valgrind (Debug)
    ccscs*)
      # Coverage (Debug) & Valgrind (Debug)
      pr=`echo $prline |  sed -r 's/^[^0-9]*([0-9]+).*/\1/'`
      logfile=$regdir/logs/ccscs-draco-Debug-coverage-master-pr${pr}.log
      allow_file_to_age $logfile 600
      echo "- Starting regression (coverage) for pr${pr}."
      echo "  Log: $logfile"
      $regdir/draco/regression/regression-master.sh -r -b Debug -e coverage \
        -p draco -f pr${pr} &> $logfile &
      logfile=$regdir/logs/ccscs-draco-Debug-valgrind-master-pr${pr}.log
      allow_file_to_age $logfile 600
      echo "- Starting regression (valgrind) for pr${pr}."
      echo "  Log: $logfile"
      $regdir/draco/regression/regression-master.sh -r -b Debug -e valgrind \
        -p draco -f pr${pr} &> $logfile &
      ;;

    # Moonlight: Fulldiagnostics (Debug)
    ml-fey*)
      pr=`echo $prline |  sed -r 's/^[^0-9]*([0-9]+).*/\1/'`
      logfile=$regdir/logs/ml-draco-Debug-fulldiagnostics-master-pr${pr}.log
      allow_file_to_age $logfile 600
      echo "- Starting regression (fulldiagnostics) for pr${pr}."
      echo "  Log: $logfile"
      $regdir/draco/regression/regression-master.sh -r -b Debug \
        -e fulldiagnostics -p draco -f pr${pr} &> $logfile &
      ;;

    # Snow ----------------------------------------
    sn-fe*)
      pr=`echo $prline |  sed -r 's/^[^0-9]*([0-9]+).*/\1/'`
      logfile=$regdir/logs/sn-draco-Debug-master-pr${pr}.log
      allow_file_to_age $logfile 600
      echo "- Starting regression for pr${pr}."
      echo "  Log: $logfile"
      $regdir/draco/regression/regression-master.sh -r -b Debug \
        -p draco -f pr${pr} &> $logfile &
      ;;

    # Trinitite: Release
    tt-fey*)
      pr=`echo $prline |  sed -r 's/^[^0-9]*([0-9]+).*/\1/'`
      logfile=$regdir/logs/tt-draco-Release-master-pr${pr}.log
      allow_file_to_age $logfile 600
      echo "- Starting regression (Release) for pr${pr}."
      echo "  Log: $logfile"
      $regdir/draco/regression/regression-master.sh -r -b Release -p draco \
        -f pr${pr} &> $logfile &
      ;;

    # Darwin ----------------------------------------
    darwin-fe*)
      # No CI
      ;;
  esac
done

# Prepare for Jayenne and Capsaicin Prs --------------------------------------
# Do we need to build draco? Only build draco-develop once per day.
eval "$(date +'today=%F now=%s')"
midnight=$(date -d "$today 0" +%s)
case ${target} in
  ccscs*) draco_tag_file=$regdir/logs/last-draco-develop-ccscs.log ;;
  ml-fe*) draco_tag_file=$regdir/logs/last-draco-develop-ml.log ;;
  sn-fe*) draco_tag_file=$regdir/logs/last-draco-develop-sn.log ;;
  tt-fe*) draco_tag_file=$regdir/logs/last-draco-develop-tt.log ;;
  darwin-fe*) draco_tag_file=$regdir/logs/last-draco-develop-darwin.log ;;
esac
draco_last_built=$(date +%s -r $draco_tag_file)

# Get the list of new Jayenne and Capsaicin Prs
jayenne_prs=`grep 'refs/merge-requests/[0-9]*/head$' $TMPFILE_JAYENNE | awk '{print $NF}'`
capsaicin_prs=`grep 'refs/merge-requests/[0-9]*/head$' $TMPFILE_CAPSAICIN | awk '{print $NF}'`

ipr=0 # count the number of PRs processed. Only the first needs to build draco.
for prline in $jayenne_prs $capsaicin_prs; do

#  seconds_since_draco_built=`expr $(date +%s) - $(date +%s -r $draco_tag_file)`

  # ----------------------------------------
  # Build draco-develop once per day for each case.
  #
  # If we haven't built draco today, build it with this PR, otherwise link to
  # the existing draco build.  Additionally, if two PRs are started at the same
  # time, only build draco for the 1st one.
  if [[ $midnight -gt $draco_last_built ]] && [[ $ipr == 0 ]]; then

    echo " "
    echo "Found a Jayenne or Capsaicin PR, but we need to build draco-develop first..."
    echo " "

    projects="draco"
    featurebranches="develop"

    # Reset the modified date on the file used to determine when draco was last
    # built.
    date &> $draco_tag_file

    case ${target} in

      # CCS-NET: Coverage (Debug) & Valgrind (Debug)
      ccscs*)
        logfile=$regdir/logs/ccscs-draco-Debug-coverage-master-develop.log
        allow_file_to_age $logfile 600
        echo "- Starting regression (coverage) for develop."
        echo "  Log: $logfile"
        $regdir/draco/regression/regression-master.sh -r -b Debug -e coverage \
          -p "${projects}" &> $logfile &

        logfile=$regdir/logs/ccscs-draco-Debug-valgrind-master-develop.log
        allow_file_to_age $logfile 600
        echo "- Starting regression (valgrind) for develop."
        echo "  Log: $logfile"
        $regdir/draco/regression/regression-master.sh -r -b Debug -e valgrind \
          -p "${projects}" &> $logfile
        # Do not put the above command into the background! It must finish
        # before jayenne is started.
        ;;

      # Moonlight: Fulldiagnostics (Debug)
      ml-fey*)
        logfile=$regdir/logs/ml-draco-Debug-fulldiagnostics-master-develop.log
        allow_file_to_age $logfile 600
        echo "- Starting regression (fulldiagnostics) for develop."
        echo "  Log: $logfile"
        $regdir/draco/regression/regression-master.sh -r -b Debug \
          -e fulldiagnostics -p draco &> $logfile
        # Do not put the above command into the background! It must finish
        # before jayenne is started.
        ;;

      # Snow ----------------------------------------
      sn-fe*)
        logfile=$regdir/logs/sn-draco-Debug-master-develop.log
        allow_file_to_age $logfile 600
        echo "- Starting regression for develop."
        echo "  Log: $logfile"
        $regdir/draco/regression/regression-master.sh -r -b Debug \
          -p draco &> $logfile
        # Do not put the above command into the background! It must finish
        # before jayenne is started.
        ;;

      # Trinitite: Release
      tt-fey*)
        logfile=$regdir/logs/tt-draco-Release-master-develop.log
        allow_file_to_age $logfile 600
        echo "- Starting regression (Release) for develop."
        echo "  Log: $logfile"
        $regdir/draco/regression/regression-master.sh -r -b Release -p draco \
          &> $logfile
        # Do not put the above command into the background! It must finish
        # before jayenne is started.
        ;;

      # Darwin ----------------------------------------
      darwin-fe*)
        # No CI
        ;;
    esac
  fi
  ((ipr++))
done

# Jayenne CI ------------------------------------------------------------

projects="jayenne"
for prline in $jayenne_prs; do

  # ----------------------------------------
  # Build Jayenne PRs against draco-develop
  #
  # All of these can be put into the backround when they run since they are
  # completely independent.
  pr=`echo $prline |  sed -r 's/^[^0-9]*([0-9]+).*/\1/'`
  featurebranches="pr${pr}"

  case ${target} in

    # CCS-NET: Coverage (Debug) & Valgrind (Debug)
    ccscs*)
      logfile=$regdir/logs/ccscs-jayenne-Debug-coverage-master-pr${pr}.log
      allow_file_to_age $logfile 600
      echo "- Starting regression (coverage) for pr${pr}."
      echo "  Log: $logfile"
      $regdir/draco/regression/regression-master.sh -r -b Debug -e coverage \
        -p "${projects}" -f "${featurebranches}" &> $logfile &

      logfile=$regdir/logs/ccscs-jayenne-Debug-valgrind-master-pr${pr}.log
      allow_file_to_age $logfile 600
      echo "- Starting regression (valgrind) for pr${pr}."
      echo "  Log: $logfile"
      $regdir/draco/regression/regression-master.sh -r -b Debug -e valgrind \
        -p "${projects}" -f "${featurebranches}" &> $logfile &
      ;;

    # Moonlight: Fulldiagnostics (Debug)
    ml-fey*)
      logfile=$regdir/logs/ml-jayenne-Debug-fulldiagnostics-master-pr${pr}.log
      allow_file_to_age $logfile 600
      echo "- Starting regression (fulldiagnostics) for pr${pr}."
      echo "  Log: $logfile"
      $regdir/draco/regression/regression-master.sh -r -b Debug \
        -e fulldiagnostics -p "${projects}" -f "${featurebranches}" \
        &> $logfile &
      ;;

    # Snow ----------------------------------------
    sn-fe*)
      logfile=$regdir/logs/sn-jayenne-Debug-master-pr${pr}.log
      allow_file_to_age $logfile 600
      echo "- Starting regression for pr${pr}."
      echo "  Log: $logfile"
      $regdir/draco/regression/regression-master.sh -r -b Debug \
        -p "${projects}" -f "${featurebranches}" \
        &> $logfile &
      ;;

    # Trinitite: Release
    tt-fey*)
      logfile=$regdir/logs/tt-jayenne-Release-master-pr${pr}.log
      allow_file_to_age $logfile 600
      echo "- Starting regression (Release) for pr${pr}."
      echo "  Log: $logfile"
      $regdir/draco/regression/regression-master.sh -r -b Release \
        -p "${projects}" -f "${featurebranches}" &> $logfile &
      ;;

    # Darwin ----------------------------------------
    darwin-fe*)
      # No CI
      ;;
  esac
done

# Capsaicin CI ------------------------------------------------------------

projects="capsaicin"
for prline in $capsaicin_prs; do

  # ----------------------------------------
  # Build Capsaicin PRs against draco-develop
  #
  # All of these can be put into the backround when they run since they are
  # completely independent.
  pr=`echo $prline |  sed -r 's/^[^0-9]*([0-9]+).*/\1/'`
  featurebranches="pr${pr}"

  case ${target} in

    # CCS-NET: Coverage (Debug) & Valgrind (Debug)
    ccscs*)
      logfile=$regdir/logs/ccscs-capsaicin-Debug-coverage-master-pr${pr}.log
      allow_file_to_age $logfile 600
      echo "- Starting regression (coverage) for pr${pr}."
      echo "  Log: $logfile"
      $regdir/draco/regression/regression-master.sh -r -b Debug -e coverage \
        -p "${projects}" -f "${featurebranches}" &> $logfile &

      logfile=$regdir/logs/ccscs-capsaicin-Debug-valgrind-master-pr${pr}.log
      allow_file_to_age $logfile 600
      echo "- Starting regression (valgrind) for pr${pr}."
      echo "  Log: $logfile"
      $regdir/draco/regression/regression-master.sh -r -b Debug -e valgrind \
        -p "${projects}" -f "${featurebranches}" &> $logfile &
      ;;

    # Moonlight: Fulldiagnostics (Debug)
    ml-fey*)
      logfile=$regdir/logs/ml-capsaicin-Debug-fulldiagnostics-master-pr${pr}.log
      allow_file_to_age $logfile 600
      echo "- Starting regression (fulldiagnostics) for pr${pr}."
      echo "  Log: $logfile"
      $regdir/draco/regression/regression-master.sh -r -b Debug \
        -e fulldiagnostics -p "${projects}" -f "${featurebranches}" \
        &> $logfile &
      ;;

    # Snow ----------------------------------------
    sn-fe*)
      logfile=$regdir/logs/sn-capsaicin-Debug-master-pr${pr}.log
      allow_file_to_age $logfile 600
      echo "- Starting regression for pr${pr}."
      echo "  Log: $logfile"
      $regdir/draco/regression/regression-master.sh -r -b Debug \
        -p "${projects}" -f "${featurebranches}" \
        &> $logfile &
      ;;

    # Trinitite: Release
    tt-fey*)
      logfile=$regdir/logs/tt-capsaicin-Release-master-pr${pr}.log
      allow_file_to_age $logfile 600
      echo "- Starting regression (Release) for pr${pr}."
      echo "  Log: $logfile"
      $regdir/draco/regression/regression-master.sh -r -b Release \
        -p "${projects}" -f "${featurebranches}" &> $logfile &
      ;;

    # Darwin ----------------------------------------
    darwin-fe*)
      # No CI
      ;;
  esac
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
run "rm $TMPFILE_DRACO $TMPFILE_JAYENNE $TMPFILE_CAPSAICIN"
run "rm $lockfile"

echo " "
echo "All done."


#------------------------------------------------------------------------------#
# End sync_repository.sh
#------------------------------------------------------------------------------#
