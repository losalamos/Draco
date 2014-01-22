#!/bin/bash -l

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

# dry_run=1

# Helpful functions:
die () { echo "ERROR: $1"; exit 1;}

run () {
   echo $1
   if ! test $dry_run; then
      eval $1
   fi
}

unpack_repo() {
   pkg=$1
   echo "Remove old files..."
   if test -f ${pkg}.hotcopy.tar; then
      run "rm -f $[pkg}.hotcopy.tar"
   fi
   if test -d ${pkg}.hotcopy; then
      run "rm -rf ${pkg}.hotcopy"
   fi
   
   echo "Unpacking SVN repository for $pkg ..."
   run "pull ${pkg}.repo ."
   run "tar -xvf ${pkg}.hotcopy.tar"
   run "mv ${pkg} ${pkg}.old"
   run "mv ${pkg}.hotcopy ${pkg}"
   run "chgrp -R draco ${pkg} ${pkg}.hotcopy.tar"
   run "chmod -R g+rwX,o-rwX ${pkg} ${pkg}.hotcopy.tar"
   echo " "
}

# working directory
start_dir=`pwd`
work_dir=/usr/projects/draco/svn
if test -d $work_dir; then
   cd $work_dir
else
   die "could not cd to $work_dir"
fi

# Ask Mercury if there are any items available for pulling from Yellow
possible_items_to_pull=`status | awk '{print $1}'`

# Loop over all items that mercury listed, if 'draco.repo' is found, then mark it for unpacking.

draco_ready=0
capsaicin_ready=0
jayenne_ready=0

for item in $possible_items_to_pull; do
   if test ${item} = "draco.repo";     then draco_ready=1; fi
   if test ${item} = "jayenne.repo";   then jayenne_ready=1; fi
   if test ${item} = "capsaicin.repo"; then capsaicin_ready=1; fi
done

# If found, pull the files
if test ${draco_ready} = 1; then unpack_repo "draco"; fi
if test ${jayenne_ready} = 1; then unpack_repo "jayenne"; fi
if test ${capsaicin_ready} = 1; then unpack_repo "capsaicin"; fi
