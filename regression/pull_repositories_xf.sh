#!/bin/bash

#------------------------------------------------------------------------------#
# Pull SVN repositories from Yellow 
#
# Assumptions:
# 1. Mercury requests must have the id tag set
#    tag             tar file              directory name
#    ---             ---------             ---------------
#    draco.repo      draco.hotcopy.tar     draco.hotcopy
#    jayenne.repo    jayenne.hotcopy.tar   jayenne.hotcopy
#    capsaicin.repo  capsaicin.hotcopy.tar capsaicin.hotcopy
# 2. SVN repositories live at /usr/projects/draco/svn
# 3. Kerberos keytab files is at $HOME/.ssh/xfkeytab and is signed
#    with principal transfer/${USER}push@lanl.gov
#------------------------------------------------------------------------------#
# How to generate keytab files:
# - See notes in push_repositories_xf.sh or the Draco Wiki.
# Obtain credentials for pull:
#   kinit -f -l 8h -kt $HOME/.ssh/xfkeytab transfer/${USER}push@lanl.gov
#------------------------------------------------------------------------------#

# dry_run=1

# Helpful functions:
die () { echo "FATAL ERROR: $1"; exit 1;}

run () {
   echo $1
   if ! test $dry_run; then
      eval $1
   fi
}

function xfpull()
{
    wantfile=$1
    filesavailable=`ssh red@transfer.lanl.gov myfiles`
    # sanity check: is the requested file in the list?
    fileready=`echo $filesavailable | grep $wantfile`
    if test "${fileready}x" = "x"; then
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

# working directory
start_dir=`pwd`
work_dir=/usr/projects/draco/svn
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

draco_ready=0
capsaicin_ready=0
jayenne_ready=0

for item in $possible_items_to_pull; do
   if test ${item} = "draco.hotcopy.tar";     then draco_ready=1; fi
   if test ${item} = "jayenne.hotcopy.tar";   then jayenne_ready=1; fi
   if test ${item} = "capsaicin.hotcopy.tar"; then capsaicin_ready=1; fi
done

# If found, pull the files
if test ${draco_ready} = 1; then unpack_repo "draco"; fi
if test ${jayenne_ready} = 1; then unpack_repo "jayenne"; fi
if test ${capsaicin_ready} = 1; then unpack_repo "capsaicin"; fi

# Update permisssions as needed
run "cd ${work_dir}/.."
run "chgrp -R draco svn"
run "chmod -R g+rwX,o-rwX svn"

# Update Module directories
cd /usr/projects/draco/vendors/environment
/usr/projects/draco/vendors/subversion-1.8.5/ml/bin/svn up

