#!/bin/bash

# Ensure that the permissions are correct
umask 0002
set MYHOSTNAME = `uname -n`

# Establish modules environment and load proper version of svn.
if test -z "$MODULESHOME"; then
   if test -f /ccs/codes/radtran/vendors/modules-3.2.9/init/bash; then
      source /ccs/codes/radtran/vendors/modules-3.2.9/init/bash
    fi
fi

# update the scripts directories in /home/regress
run () {
    echo $1
    if ! [ $dry_run ]; then
       eval $1
    fi
}

run "module use /ccs/codes/radtran/vendors/Modules"
run "module load svn"

run "cd /home/regress/draco/config; svn update"
run "cd /home/regress/draco/regression; svn update"
run "cd /home/regress/draco/environment; svn update"
run "cd /home/regress/jayenne/regression; svn update"
run "cd /home/regress/capsaicin/scripts; svn update"
run "cd /home/regress/asterisk/regression; svn update"

