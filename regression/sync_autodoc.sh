#!/bin/bash
#!/bin/bash
##---------------------------------------------------------------------------##
## File  : sync_autodoc.sh
## Date  : Wednesday, Oct 17, 2018, 12:06 pm
## Author: Kelly Thompson <kgt@lanl.gov>
## Note  : Copyright (C) 2018, Los Alamos National Security, LLC.
##         All rights are reserved.
##
## Look for files at /ccs/codes/radtran/autodoc and rsync (publish) them to the
## https://rtt.lanl.gov/autodoc web server.
##---------------------------------------------------------------------------##

umask 0002
function run ()
{
    echo $1;
    if ! [ $dry_run ]; then
        eval $1;
    fi
}

# Sanity checks:
mach=`uname -n`
if test "${mach}" != "ccsnet3.lanl.gov"; then
   echo "FATAL ERROR: This script must be run from rtt.lanl.gov."
   exit 1
fi

# Locations
sourcedir=/ccs/codes/radtran/autodoc
destdir=/var/www/virtual/rtt.lanl.gov/html/autodoc

if ! test -d $sourcedir; then
   echo "FATAL ERROR: Cannot find sourcedir = $sourcedir"
   exit 1
fi
if ! test -d $destdir; then
   echo "FATAL ERROR: Cannot find destdir = $destdir"
   exit 1
fi

# Remove old files and directories.
current_dir=`pwd`

echo -e "\nCleaning up old builds."
run "cd $sourcedir"
# run "find . -mtime +5 -type f"
# run "find . -mtime +5 -type f -delete"
run "find . -maxdepth 1 -mtime +3 -name 'pr*' -type d"
run "find . -maxdepth 1 -mtime +3 -name 'pr*' -type d -exec rm -rf {} \;"
echo " "
run "cd $destdir"
# run "find . -mtime +5 -type f"
# run "find . -mtime +5 -type f -delete"
run "find . -maxdepth 1 -mtime +3 -name 'pr*' -type d"
run "find . -maxdepth 1 -mtime +3 -name 'pr*' -type d -exec rm -rf {} \;"

# Copy new files and directory to the published location
run "cd $current_dir"
run "rsync --delete -vaum ${sourcedir}/ ${destdir}"
run "chgrp -R draco ${destdir}"
run "chmod -R g+rwX,o=g-w ${destdir}"

#------------------------------------------------------------------------------#
# End sync_autodoc.sh
#------------------------------------------------------------------------------#
