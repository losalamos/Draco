#!/bin/tcsh -f

# This script is used to mirror portions of the Draco/Jayenne SVN
# repositories to a HPC location.  The repository must be mirrored
# because the ctest regression system must be run from the HPC backend
# (via msub) where access to ccscs8:/ccs/codes/radtran/svn is not
# available.

# source system dotfiles if running batch (cron) and not already set up
# if( -f /etc/csh.cshrc ) then
#   # source /etc/csh.cshrc
#   if( -f /etc/csh.login ) then
#     source /etc/csh.login
#   endif
# endif

# # source file to set up module alias if not already set
# set _is_set=`alias module`
# if ( "$_is_set" == "" ) then
#   foreach _loc ( /opt/modules/default/etc /etc/profile.d )
#     if ( -e $_loc/modules.csh ) then
#       source $_loc/modules.csh
#       break
#     endif
#   end
# endif

# Ensure that the permissions are correct
umask 0002
set MYHOSTNAME = `uname -n`

# Credentials via Keychain (SSH)
# http://www.cyberciti.biz/faq/ssh-passwordless-login-with-keychain-for-scripts
/usr/projects/draco/vendors/keychain-2.7.1/keychain $HOME/.ssh/cmake_dsa
if ( -f $HOME/.keychain/$MYHOSTNAME.-csh ) then
    source $HOME/.keychain/$MYHOSTNAME.-csh
endif
if ( -f $HOME/.keychain/$MYHOSTNAME-csh ) then
    source $HOME/.keychain/$MYHOSTNAME-csh
endif

# # Ensure correct svn is available in the environment
# if ( -d /users/kellyt/draco/environment/Modules ) then
#     module use /users/kellyt/draco/environment/Modules/hpc
#     module use /users/kellyt/draco/environment/Modules/ct-fe
#     module load svn
# endif
# if ( "$MYHOSTNAME" == "ct-fe1" ) then
#     module use /usr/projects/hpcsoft/modulefiles/cielito/hpc-tools
#     module load subversion
# endif


# update the scripts directories in /usr/projects/jayenne
echo "cd /usr/projects/jayenne/regress/draco/config; svn update"
cd /usr/projects/jayenne/regress/draco/config
svn update
echo "cd /usr/projects/jayenne/regress/draco/regression; svn update"
cd /usr/projects/jayenne/regress/draco/regression
svn update
echo "cd /usr/projects/jayenne/regress/Modules; svn update"
cd /usr/projects/jayenne/regress/Modules
svn update

# SVN portions
# ------------------------------------------------------------

set svnroot = /usr/projects/jayenne/regress/svn

# Setup directory structure
if ( -d $svnroot ) then
    :
else
    echo "*** ERROR ***"
    echo "*** SVN repository not found ***"
    exit 1
    # http://journal.paul.querna.org/articles/2006/09/14/using-svnsync/
    # mkdir -p ${svnroot}
    # svnadmin create ${svnroot}/jayenne
    # (update hooks/pre-commit and hooks/pre-svnprop; chmod 775 )
    # svnsync init file:///${svnroot}/jayenne svn+ssh://ccscs8/ccs/codes/radtran/svn/jayenne
    # svnsync sync file:///${svnroot}/jayenne
endif

# Draco
echo "svnsync --non-interactive sync file:///${svnroot}/draco"
svnsync --non-interactive sync file:///${svnroot}/draco

# Jayenne
echo "svnsync --non-interactive sync file:///${svnroot}/jayenne"
svnsync --non-interactive sync file:///${svnroot}/jayenne
#chgrp -R jayenne ${svnroot}

# Capsaicin
echo "svnsync --non-interactive sync file:///${svnroot}/capsaicin"
svnsync --non-interactive sync file:///${svnroot}/capsaicin

# also update the scripts directory
if ( -d /usr/projects/jayenne/regress/jayenne/regression ) then
    echo "cd /usr/projects/jayenne/regress/jayenne/regression; svn update"
    cd /usr/projects/jayenne/regress/jayenne/regression
    svn update
else
    echo "mkdir -p /usr/projects/jayenne/regress/jayenne"
    mkdir -p /usr/projects/jayenne/regress/jayenne
    echo "svn co svn+ssh://ccscs8/ccs/codes/radtran/svn/jayenne/jayenne-project/regression regression"
    svn co svn+ssh://ccscs8/ccs/codes/radtran/svn/jayenne/jayenne-project/regression regression
endif

if ( -d /usr/projects/jayenne/regress/capsaicin/scripts ) then
    echo "cd /usr/projects/jayenne/regress/capsaicin/scripts; svn update"
    cd /usr/projects/jayenne/regress/capsaicin/scripts
    svn update
else
    echo "mkdir -p /usr/projects/jayenne/capsaicin/scripts"
    mkdir -p /usr/projects/jayenne/capsaicin/scripts
    echo "svn co svn+ssh://ccscs8/ccs/codes/radtran/svn/capsaicin/trunk/scripts scripts"
    svn co svn+ssh://ccscs8/ccs/codes/radtran/svn/capsaicin/trunk/scripts scripts
endif


# Notes:
# ------------------------------------------------------------
# McKay uses the following:
# Purpose: Get Kerberos ticket
# Command: /usr/kerberos/bin/kinit -k -t /users/lmdm/.ssh/keytabfile lmdm@lanl.gov
# Output:  
