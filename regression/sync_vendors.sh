#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/sync_vendors.sh
## Date  : Tuesday, Oct 25, 2016, 09:07 am
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
lockfile=/var/tmp/sync_vendors_$target.lock
[ "${FLOCKER}" != "${lockfile}" ] && exec env FLOCKER="${lockfile}" flock -en "${lockfile}" "${0}" "$@" || :

# All output will be saved to this log file.  This is also the lockfile for flock.
timestamp=`date +%Y%m%d-%H%M`
logdir="$( cd $scriptdir/../../logs && pwd )"
logfile=$logdir/sync_vendors-$target-$timestamp.log

# Redirect all output to a log file.
exec > $logfile
exec 2>&1

# import some bash functions
source $scriptdir/scripts/common.sh

echo -e "Executing $0 $*...\n"
echo "Group: `id -gn`"
echo "umask: `umask`"

#------------------------------------------------------------------------------#
# Run from ccscs7:
#   - Mirror ccscs7:/scratch/vendors -> /ccs/codes/radtran/vendors/rhel72
#     Do not mirror the ndi directory (300+ GB)
#   - Mirror ccscs7:/scratch/vendors  -> ccscs[234568]:/scratch/vendors
#     Do not mirror ndi to ccscs5 due to limited size of /scratch.

# From dir (ccscs7)
vdir=/scratch/vendors
# To dir
# r72v=/ccs/codes/radtran/vendors/rhel72vendors
# To machines (cccs[234568]:/scratch/vendors)
# omit ccscs5 (scratch is too full)
ccs_servers="ccscs1 ccscs2 ccscs3 ccscs6 ccscs8 ccscs9"

# exclude these directories
exclude_dirs="spack.mirror spack.rasa tmp spack.temporary spack.test spack.ccs.developmental"
exclude_items=""
for item in $exclude_dirs; do
  exclude_items+=" --exclude=$item"
done

# Credentials via Keychain (SSH)
# http://www.cyberciti.biz/faq/ssh-passwordless-login-with-keychain-for-scripts
#$vdir/keychain-2.8.2/keychain $HOME/.ssh/cmake_dsa $HOME/.ssh/cmake_rsa
$vdir/keychain-2.8.2/keychain $HOME/.ssh/id_rsa
if test -f $HOME/.keychain/$HOSTNAME-sh; then
  run "source $HOME/.keychain/$HOSTNAME-sh"
else
  echo "Error: could not find $HOME/.keychain/$HOSTNAME-sh"
fi

# Banner
echo "Rsync /scratch/vendors to /scratch/vendors on:"
for m in $ccs_servers; do echo " - ${m}"; done

# Sanity check
if ! test -d $vdir; then
  echo "Source directory $vdir is missing."
  exit 1
fi
if ! test -d ${vdir}-ec; then
  echo "Source directory ${vdir}-ec is missing."
  exit 1
fi

# Reset permissions (if possible):
echo " "
echo "Clean up permissions on source files..."
run "cd $vdir"
vdir_subdirs=`\ls -1 $vdir`
for dir in $vdir_subdirs; do
  run "chgrp -R draco $dir &> /dev/null"
  run "chmod -R g+rX,o+rX $dir &> /dev/null"
  run "find $dir -type d -exec chmod g+s {} \;"
done

run "cd ${vdir}-ec"
vdir_subdirs=`\ls -1 ${vdir}-ec`
for dir in $vdir_subdirs; do
  run "chgrp -R ccsrad $dir &> /dev/null"
  run "chmod -R g+rX,o-rwX $dir &> /dev/null"
  run "find $dir -type d -exec chmod g+s {} \;"
done

#echo " "
#echo "Save a copy of /scratch/vendors to $r72v..."
#run "rsync -av --exclude 'ndi' --delete $vdir/ $r72v"
# rsync -av --omit-dir-times --checksum --human-readable --progress <local dir> <remote dir>

# rsync vendors ccscs7 -> other machines.
# but limit network to 50 MB/sec (400 mbps)
echo " "
echo "Rsync $vdir and ${vdir}-ec to other ccs-net servers... "
for m in $ccs_servers; do
  if [[ `uname -n | grep -c $m` = 0 ]]; then
    case $m in
      ccscs5)
        # do not copy ndi, not enough space
        run "rsync -av ${exclude_items} --exclude=ndi --delete --bwlimit=50000 $vdir/ ${m}:$vdir"
        run "rsync -av ${exclude_items} --exclude=ndi --delete --bwlimit=50000 ${vdir}-ec/ ${m}:${vdir}-ec"
        ;;
      *)
        # run "rsync -av --delete --bwlimit=50000 $vdir/ ${m}:$vdir"
        run "rsync -av ${exclude_items} --delete --bwlimit=200000 $vdir/ ${m}:$vdir"
        run "rsync -av ${exclude_items} --delete --bwlimit=200000 ${vdir}-ec/ ${m}:${vdir}-ec"
        ;;
    esac
  fi
done

# Cleanup
echo " "
echo "Cleaning up..."
run "rm $lockfile"

echo -e "\n--------------------------------------------------------------------------------"
echo "All done."
echo "--------------------------------------------------------------------------------"

#------------------------------------------------------------------------------#
# End sync_vendors.sh
#------------------------------------------------------------------------------#
