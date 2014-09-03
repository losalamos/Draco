#!/bin/bash

umask 0002

target="`uname -n | sed -e s/[.].*//`"
arch=`uname -m`

# Ensure that the permissions are correct
case ${target} in
darwin* | cn[0-9]*)
   SVN=/projects/opt/subversion/1.7.14/bin/svn
   REGDIR=/projects/opt/draco/regress
   ;;
*)
   SVN=/ccs/codes/radtran/vendors/subversion-1.8.5/bin/svn
   REGDIR=/home/regress
   ;;
esac

# Helper function
run () {
    echo $1
    if ! [ $dry_run ]; then eval $1; fi
}

run "cd ${REGDIR}/draco/config; ${SVN} update"
run "cd ${REGDIR}/draco/regression; ${SVN} update"
run "cd ${REGDIR}/draco/environment; ${SVN} update"
run "cd ${REGDIR}/jayenne/regression; ${SVN} update"
run "cd ${REGDIR}/capsaicin/scripts; ${SVN} update"
run "cd ${REGDIR}/asterisk/regression; ${SVN} update"

