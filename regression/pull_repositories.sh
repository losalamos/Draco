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
# 3. Kerberos keytab files is at $HOME/.ssh/cron.keytab and is signed
#    with principal $USER@lanl.gov
#------------------------------------------------------------------------------#
# How to generate keytab files:
#
# /usr/kerberos/sbin/kadmin -p ${USER}@lanl.gov
# > ktadd -k $HOME/.ssh/cron.keytab -p ${USER}@lanl.gov
# > quit
# chmod 600 $HOME/.ssh/cron.keytab

# Enable kerberos credentials via keytab authentication:
#
# kinit -f -l 8h -kt $HOME/.ssh/cron.keytab ${USER}@lanl.gov
#------------------------------------------------------------------------------#

STATUS=/usr/local/bin/status
PULL=/usr/local/bin/pull

# dry_run=1

# Helpful functions:
die () { echo "FATAL ERROR: $1"; exit 1;}

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
      run "rm -f ${pkg}.hotcopy.tar"
   fi
   if test -d ${pkg}.hotcopy; then
      run "rm -rf ${pkg}.hotcopy"
   fi
   if test -d ${pkg}.old; then
      run "rm -rf ${pkg}.old"
   fi
   
   echo "Unpacking SVN repository for $pkg ..."
   run "${PULL} ${pkg}.repo ."
   run "tar -xvf ${pkg}.hotcopy.tar"
   run "mv ${pkg} ${pkg}.old"
   run "mv ${pkg}.hotcopy ${pkg}"
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

# Ensure we have a kerberos ticket
run "kinit -f -l 8h -kt $HOME/.ssh/cron.keytab ${USER}@lanl.gov"

# Ask Mercury if there are any items available for pulling from Yellow
possible_items_to_pull=`${STATUS} | awk '{print $1}'`

# Loop over all items that mercury listed, if 'draco.repo' is found,
# then mark it for unpacking. 

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

# Update permisssions as needed
run "cd ${work_dir}/.."
run "chgrp -R draco svn"
run "chmod -R g+rwX,o-rwX svn"
