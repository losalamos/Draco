#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/sync_vendors.sh
## Date  : Tuesday, Oct 25, 2016, 09:07 am
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# switch to group 'ccsrad' and set umask
if [[ $(id -gn) != ccsrad ]]; then
  exec sg ccsrad "$0 $*"
fi
umask 0002

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
r72v=/ccs/codes/radtran/vendors/rhel72vendors
# To machines (cccs[234568]:/scratch/vendors)
# omit ccscs5 (scratch is too full)
ccs_servers="ccscs2 ccscs3 ccscs4 ccscs6 ccscs8"

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
run "chgrp -R draco $vdir &> /dev/null"
run "chmod -R g+rwX,o=g-w $vdir &> /dev/null"
run "chmod -R o-rwX $vdir/ndi* $vdir/csk* $vdir/cubit* $vdir/eospac* &> /dev/null"

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

echo -e "\n--------------------------------------------------------------------------------"
echo "All done."
echo "--------------------------------------------------------------------------------"

#------------------------------------------------------------------------------#
# End sync_vendors.sh
#------------------------------------------------------------------------------#
