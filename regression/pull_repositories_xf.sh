#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/pull_repositories_xf.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2018, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##
# Pull Git repositories from Yellow
#
# Assumptions:
# 1. Tar'ed repository names (listed at repository_lsit.sh)
#    repo             tar file
#    ---             ---------
#    Draco           lanl_Draco.git.tar
#    branson         lanl_branson.git.tar
#    jayenne         jayenne_jayenne.git.tar
#    imcdoc          jayenne_imcdoc.git.tar
#    capsaicin       capsaicin_capsaicin.git.tar
#    capsaicin/core  capsaicin_core.git.tar
#    capsaicin/docs  capsaicin_docs.git.tar
#    capsaicin/npt   capsaicin_npt.git.tar
#    capsaicin/trt   capsaicin_trt.git.tar
# 2. Git repositories live at /usr/projects/draco/git
# 3. Kerberos keytab files is at $HOME/.ssh/xfkeytab and is signed
#    with principal transfer/${USER}push@lanl.gov
#------------------------------------------------------------------------------#
# How to generate keytab files:
# - See notes in push_repositories_xf.sh or the Draco Wiki.
# Obtain credentials for pull:
#   kinit -f -l 8h -kt $HOME/.ssh/xfkeytab transfer/${USER}push@lanl.gov
#------------------------------------------------------------------------------#

# dry_run=1

# load some common bash functions
export scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [[ -f $scriptdir/scripts/common.sh ]]; then
  source $scriptdir/scripts/common.sh
else
  echo " "
  echo "FATAL ERROR: Unable to locate Draco's bash functions: "
  echo "   looking for .../regression/scripts/common.sh"
  echo "   searched scriptdir = $scriptdir"
  exit 1
fi

#------------------------------------------------------------------------------#
function xfpull()
{
  wantfile=$1
  filesavailable=`ssh red@transfer.lanl.gov myfiles`
  # sanity check: is the requested file in the list?
  fileready=`echo $filesavailable | grep $wantfile`
  if [[ ! ${fileready} ]]; then
    echo "ERROR: File '${wantfile}' is not available (yet?) to pull."
    echo "       Run 'xfstatus' to see list of available files."
    return
  fi
  # Find the file identifier for the requested file.  The variable
  # filesavailable contains a list of pairs:
  # { (id1, file1), (id2, file2), ... }.  Load each pair and if the
  # filename matches the requested filename then pull that file id.
  # Once pulled, remove the file from transfer.lanl.gov.
  is_file_id=1
  for entry in $filesavailable; do
    if test $is_file_id = 1; then
      fileid=$entry
      is_file_id=0
    else
      if test $entry = $wantfile; then
        echo "scp red@transfer.lanl.gov:${fileid} ."
        scp red@transfer.lanl.gov:${fileid} .
        echo "ssh red@transfer.lanl.gov delete ${fileid}"
        ssh red@transfer.lanl.gov delete ${fileid}
        return
      fi
      is_file_id=1
    fi
  done
}

#------------------------------------------------------------------------------#
unpack_repo_git() {
  if [[ $# != 2 ]]; then
    die "wrong number of parameters given to 'unpack_repo_git'. Expecting two arguments".
  fi
  namespace=$1
  repo=$2
  # pkg should have the form ${namespace}_${repo}.git.tar

  echo -e "\nRemove old files/directories...\n"
  if [[ -f ${namespace}_${repo}.git.tar ]]; then
    run "rm -f ${namespace}_${repo}.git.tar"
  fi
  if [[ -d ${namespace}/${repo}.git ]]; then
    if [[ -d ${namespace}/${repo}.git.old ]]; then
      run "rm -rf ${namespace}/${repo}.git.old"
    fi
  fi

  echo -e "\nUnpacking GIT repository for ${namespace}/${repo}.git...\n"
  run "xfpull ${namespace}_${repo}.git.tar"
  if [[ -d ${namespace}/${repo}.git ]]; then
    run "mv ${namespace}/${repo}.git ${namespace}/${repo}.git.old"
   fi
  run "tar -xvf ${namespace}_${repo}.git.tar"
}

#------------------------------------------------------------------------------#
# working directory
start_dir=`pwd`
gitroot=/usr/projects/draco/git
if test -d $gitroot; then
   run "cd $gitroot"
else
   die "could not cd to $gitroot"
fi

# Ensure we have a kerberos ticket
run "kinit -f -l 1h -kt $HOME/.ssh/xfkeytab transfer/${USER}push@lanl.gov"

# Ask Mercury if there are any items available for pulling from Yellow
possible_items_to_pull=`ssh red@transfer.lanl.gov myfiles | awk '{print $2}'`

# List of repositories (also used by sync_repositories.sh and
# pull_repositories_xf.sh).  It defines $git_projects.
source ${scriptdir}/repository_list.sh

# Loop over each known repository.  If a new tar file is found in the transfer
# system, pull it and unpack it at $gitroot.
for project in ${git_projects[@]}; do

  echo -e "\nProcessing $project...\n"

  namespace=`echo $project | sed -e 's%/.*%%'`
  repo=`echo $project | sed -e 's%.*/%%'`
  xf_file=${namespace}_${repo}.git.tar

  # keep track of how many files transfer knows about.  Sometimes there are
  # extra files because a machine was down when the cronjob normally runs.
  num_avail=`echo ${possible_items_to_pull} | grep -c ${xf_file}`
  while [[ $num_avail -gt 0 ]]; do
    let "num_avail -= 1"
    unpack_repo_git "${namespace}" "${repo}"
    # sometimes it takes some time before the file is removed from transfer.
    # So, if we need to unpack the tar file again, then pause first to ensure we
    # are unpacking the 2nd file in the list instead of repeating what we just
    # did.
    if [[ $num_avail -gt 0 ]]; then sleep 30; fi
  done

done

# Update permisssions as needed
run "cd ${gitroot}/.."
run "chgrp -R draco git"
run "chmod -R g+rwX,o-rwX git"
run "cd $start_dir"

echo -e "\nAll done.\n"

#------------------------------------------------------------------------------#
# End pull_repositories_xf.sh
#------------------------------------------------------------------------------#
