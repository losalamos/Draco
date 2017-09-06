#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/pull_repositories_xf.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##
# Pull Git repositories from Yellow
#
# Assumptions:
# 1. Tar'ed repository names
#    repo             tar file
#    ---             ---------
#    draco           draco.git.tar
#    jayenne         jayenne.git.tar
#    capsaicin       capsaicin.git.tar
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
export rscriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [[ -f $rscriptdir/scripts/common.sh ]]; then
  source $rscriptdir/scripts/common.sh
else
  echo " "
  echo "FATAL ERROR: Unable to locate Draco's bash functions: "
  echo "   looking for .../regression/scripts/common.sh"
  echo "   searched rscriptdir = $rscriptdir"
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
unpack_repo() {
   pkg=$1
   echo "Remove old files..."
   if test -f ${pkg}.hotcopy.tar; then
      run "rm -f ${pkg}.hotcopy.tar"
   fi
   if test -d ${pkg}.hotcopy; then
      run "rm -rf ${pkg}.hotcopy"
   fi
   if test -d ${pkg}; then
     if test -d ${pkg}.old; then
        run "rm -rf ${pkg}.old"
     fi
   fi

   echo "Unpacking SVN repository for $pkg ..."
   run "xfpull ${pkg}.hotcopy.tar"
   run "tar -xvf ${pkg}.hotcopy.tar"
   if test -d ${pkg}.hotcopy; then
      run "mv ${pkg} ${pkg}.old"
      run "mv ${pkg}.hotcopy ${pkg}"
   fi
   echo " "
}

#------------------------------------------------------------------------------#
unpack_repo_git() {
  pkg=$1
  echo "Remove old files..."
  if test -f ${pkg}.tar; then
    run "rm -f ${pkg}.tar"
  fi
  if test -d ${pkg}; then
    if test -d ${pkg}.old; then
      run "rm -rf ${pkg}.old"
    fi
  fi

  echo "Unpacking GIT repository for $pkg ..."
  run "xfpull ${pkg}.tar"
  if test -d ${pkg}; then
    run "mv ${pkg} ${pkg}.old"
   fi
  run "tar -xvf ${pkg}.tar"
  echo " "
}

#------------------------------------------------------------------------------#
# working directory
start_dir=`pwd`
work_dir=/usr/projects/draco/git
if test -d $work_dir; then
   run "cd $work_dir"
else
   die "could not cd to $work_dir"
fi

# Ensure we have a kerberos ticket
run "kinit -f -l 1h -kt $HOME/.ssh/xfkeytab transfer/${USER}push@lanl.gov"

# Ask Mercury if there are any items available for pulling from Yellow
possible_items_to_pull=`ssh red@transfer.lanl.gov myfiles | awk '{print $2}'`

# Loop over all items that mercury listed, if 'draco.repo' is found,
# then mark it for unpacking.

capsaicin_ready=0
jayenne_ready=0
draco_git_ready=0

for item in $possible_items_to_pull; do
   if test ${item} = "capsaicin.git.tar"; then capsaicin_ready=1; fi
   if test ${item} = "Draco.git.tar";     then draco_git_ready=1; fi
   if test ${item} = "jayenne.git.tar";   then jayenne_ready=1;   fi
done

# If found, pull the files
run "cd ${work_dir}"
if test ${draco_git_ready} = 1; then unpack_repo_git "Draco.git"; fi
if test ${jayenne_ready} = 1; then unpack_repo_git "jayenne.git"; fi
if test ${capsaicin_ready} = 1; then unpack_repo_git "capsaicin.git"; fi

# Update permisssions as needed
run "cd ${work_dir}/.."
run "chgrp -R draco git"
run "chmod -R g+rwX,o-rwX git"

#------------------------------------------------------------------------------#
# End pull_repositories_xf.sh
#------------------------------------------------------------------------------#
