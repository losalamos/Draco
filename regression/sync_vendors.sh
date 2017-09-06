#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/sync_vendors.sh
## Date  : Tuesday, Oct 25, 2016, 09:07 am
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# Run from ccscs7:
#   - Mirror ccscs7:/scratch/vendors -> /ccs/codes/radtran/vendors/rhel72
#     Do not mirror the ndi directory (300+ GB)
#   - Mirror ccscs7:/scratch/vendors  -> ccscs[234568]:/scratch/vendors
#     Do not mirror ndi to ccscs5 due to limited size of /scratch.

# From dir (ccscs7)
vdir=/scratch/vendors
# To dir
r72v=/ccs/codes/radtran/vendors/rhel72vendors
# To machines (cccs[234568]:/scratch/vendors)
ccs_servers="ccscs2 ccscs3 ccscs4 ccscs5 ccscs6 ccscs8"

target="`uname -n | sed -e s/[.].*//`"

# Locate the directory that this script is located in:
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# All output will be saved to this log file.  This is also the lockfile for flock.
logdir="$( cd $scriptdir/../../logs && pwd )"
logfile=$logdir/sync_vendors_$target.log
lockfile=/var/tmp/sync_vendors_$target.lock

# Prevent multiple copies of this script from running at the same time:
[ "${FLOCKER}" != "${lockfile}" ] && exec env FLOCKER="${lockfile}" flock -en "${lockfile}" "${0}" "$@" || :

# Redirect all output to a log file.
exec > $logfile
exec 2>&1

# import some bash functions
source $scriptdir/scripts/common.sh

# Ensure that the permissions are correct
run "umask 0002"

# Credentials via Keychain (SSH)
# http://www.cyberciti.biz/faq/ssh-passwordless-login-with-keychain-for-scripts
$vdir/keychain-2.8.2/keychain $HOME/.ssh/cmake_dsa $HOME/.ssh/cmake_rsa
if test -f $HOME/.keychain/$HOSTNAME-sh; then
  run "source $HOME/.keychain/$HOSTNAME-sh"
else
  echo "Error: could not find $HOME/.keychain/$HOSTNAME-sh"
fi

# Banner
echo "Rsync /scratch/vendors to $r72v and to /scratch/vendors on:"
for m in $ccs_servers; do echo " - ${m}"; done

# Sanity check
if ! test -d $vdir; then
  echo "Source directory $vdir is missing."
  exit 1
fi

# Make a backup copy of vendors to $r72v
echo " "
echo "Clean up permissions on source files..."
run "chgrp -R draco $vdir"
run "chmod -R g+rwX,o=g-w $vdir"
run "chmod -R o-rwX $vdir/ndi* $vdir/csk* $vdir/cubit* $vdir/eospac*"

echo " "
echo "Save a copy of /scratch/vendors to $r72v..."
#run "rsync -av --exclude 'ndi' --delete $vdir/ $r72v"
# rsync -av --omit-dir-times --checksum --human-readable --progress <local dir> <remote dir>

# rsync vendors ccscs7 -> other machines.
# but limit network to 50 MB/sec (400 mbps)
echo " "
echo "Rsync $vdir to other ccs-net servers... "
for m in $ccs_servers; do
  if [[ `uname -n | grep -c $m` = 0 ]]; then
    case $m in
      ccscs5)
        # do not copy ndi, not enough space
        run "rsync -av --exclude ndi --delete --bwlimit=50000 $vdir/ ${m}:$vdir"
        ;;
      *)
        run "rsync -av --delete --bwlimit=50000 $vdir/ ${m}:$vdir"
        ;;
    esac
  fi
done

# Cleanup
echo " "
echo "Cleaning up..."
run "rm $lockfile"

echo " "
echo "All done."

#------------------------------------------------------------------------------#
# End sync_vendors.sh
#------------------------------------------------------------------------------#
