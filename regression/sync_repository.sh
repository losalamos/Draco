#!/bin/bash

# This script is used to mirror portions of the Draco/Jayenne SVN
# repositories to a HPC location.  The repository must be mirrored
# because the ctest regression system must be run from the HPC backend
# (via msub) where access to ccscs8:/ccs/codes/radtran/svn is not
# available.

#
# MODULES
#
# Determine if the module command is available
modcmd=`declare -f module`
# If not found, look for it in /usr/share/Modules (ML)
if test "${modcmd}x" = "x"; then
   if test -f /usr/share/Modules/init/bash; then
      source /usr/share/Modules/init/bash
   else
      echo "ERROR: The module command was not found. No modules will be loaded."
   fi
fi
modcmd=`declare -f module`

#
# Environment
#

run () {
    echo $1
    if ! test $dry_run; then eval $1; fi
}

run "module load user_contrib svn"

# Ensure that the permissions are correct
run "umask 0002"
MYHOSTNAME=`uname -n`
regdir=/usr/projects/jayenne/regress
svnroot=$regdir/svn
svnhostmachine=ccscs7

# Credentials via Keychain (SSH)
# http://www.cyberciti.biz/faq/ssh-passwordless-login-with-keychain-for-scripts
/usr/projects/draco/vendors/keychain-2.7.1/keychain $HOME/.ssh/cmake_dsa
if test -f $HOME/.keychain/$MYHOSTNAME-sh; then
    run "source $HOME/.keychain/$MYHOSTNAME-sh"
else
    echo "Error: could not find $HOME/.keychain/$MYHOSTNAME-sh"
fi

#
# update the scripts directories in /usr/projects/jayenne
#
run "cd /usr/projects/draco/vendors/Modules; svn update"

dirs="jayenne/regression capsaicin/scripts asterisk/regression"

if test -d $regdir/draco/config; then
    run "cd $regdir/draco/config; svn update"
else
    run "mkdir -p $regdir/draco; cd $regdir/draco"
    run "svn co svn+ssh://$svnhostmachine/ccs/codes/radtran/svn/draco/trunk/config"
fi

if test -d $regdir/draco/regression; then
    run "cd $regdir/draco/regression; svn update"
else
    run "mkdir -p $regdir/draco; cd $regdir/draco"
    run "svn co svn+ssh://$svnhostmachine/ccs/codes/radtran/svn/draco/trunk/regression"
fi

if test -d $regdir/draco/environment; then
    run "cd $regdir/draco/environment; svn update"
else
    run "mkdir -p $regdir/draco; cd $regdir/draco"
    run "svn co svn+ssh://$svnhostmachine/ccs/codes/radtran/svn/draco/trunk/environment"
fi

if test -d $regdir/jayenne/regression; then
    run "cd $regdir/jayenne/regression; svn update"
else
    run "mkdir -p $regdir/jayenne; cd $regdir/jayenne"
    run "svn co svn+ssh://$svnhostmachine/ccs/codes/radtran/svn/jayenne-project/regression"
fi

if test -d $regdir/capsaicin/scripts; then
    run "cd $regdir/capsaicin/scripts; svn update"
else
    run "mkdir -p $regdir/capsaicin; cd $regdir/capsaicin"
    run "svn co svn+ssh://$svnhostmachine/ccs/codes/radtran/svn/capsaicin/trunk/scripts"
fi

if test -d $regdir/asterisk/regression; then
    run "cd $regdir/asterisk/regression; svn update"
else
    run "mkdir -p $regdir/asterisk; cd $regdir/asterisk"
    run "svn co svn+ssh://$svnhostmachine/ccs/codes/radtran/svn/asterisk/trunk/regression"
fi

#
# Sync the repository to $svnroot
#
if ! test -d $svnroot; then
    echo "*** ERROR ***"
    echo "*** SVN repository not found ***"
    exit 1
    # http://journal.paul.querna.org/articles/2006/09/14/using-svnsync/
    # mkdir -p ${svnroot}; cd ${svnroot}
    # svnadmin create ${svnroot}/jayenne
    # chgrp -R draco jayenne; chmod -R g+rwX,o=g-w jayenne
    # cd jayenne/hooks
    # cp pre-commit.tmpl pre-commit; chmod 775 pre-commit
    # vi pre-commit; comment out all code and add...
    #if ! test `whoami` = 'kellyt'; then
    #echo "This is a read only repository.  The real SVN repository is"
    #echo "at svn+ssh://ccscs8/ccs/codes/radtran/svn/draco."
    #exit 1
    #fi
    #exit 0
    # cp pre-revprop-change.tmpl pre-revprop-change; chmod 775 \
    #    pre-revprop-change
    # vi pre-revprop-change --> comment out all code.
    # cd $svnroot
    # svnsync init file:///${svnroot}/jayenne svn+ssh://ccscs8/ccs/codes/radtran/svn/jayenne
    # svnsync sync file:///${svnroot}/jayenne
fi

run "svnsync --non-interactive sync file:///${svnroot}/draco"
run "svnsync --non-interactive sync file:///${svnroot}/jayenne"
run "svnsync --non-interactive sync file:///${svnroot}/capsaicin"
run "svnsync --non-interactive sync file:///${svnroot}/asterisk"
