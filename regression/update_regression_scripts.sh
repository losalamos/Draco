#!/bin/bash

# Ensure that the permissions are correct
umask 0002
SVN=/ccs/codes/radtran/vendors/subversion-1.8.5/bin/svn

# Helper function
run () {
    echo $1
    if ! [ $dry_run ]; then
       eval $1
    fi
}

run "cd /home/regress/draco/config; ${SVN} update"
run "cd /home/regress/draco/regression; ${SVN} update"
run "cd /home/regress/draco/environment; ${SVN} update"
run "cd /home/regress/jayenne/regression; ${SVN} update"
run "cd /home/regress/capsaicin/scripts; ${SVN} update"
run "cd /home/regress/asterisk/regression; ${SVN} update"

