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
#     module use /usr/projects/jayenne/regress/draco/environment/Modules/hpc
#     module use /usr/projects/jayenne/regress/draco/environment/Modules/tu-fe
#     module load svn
# endif
# if ( "$MYHOSTNAME" == "ct-fe1" ) then
#     module use /usr/projects/hpcsoft/modulefiles/cielito/hpc-tools
#     module load subversion
# endif


# update the scripts directories in /usr/projects/jayenne
echo "cd /usr/projects/draco/vendors/Modules; svn update"
cd /usr/projects/draco/vendors/Modules
svn update
echo "cd /usr/projects/jayenne/regress/draco/config; svn update"
cd /usr/projects/jayenne/regress/draco/config
svn update
echo "cd /usr/projects/jayenne/regress/draco/regression; svn update"
cd /usr/projects/jayenne/regress/draco/regression
svn update
echo "cd /usr/projects/jayenne/regress/draco/environment; svn update"
cd /usr/projects/jayenne/regress/draco/environment
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

# Asterisk
echo "svnsync --non-interactive sync file:///${svnroot}/asterisk"
svnsync --non-interactive sync file:///${svnroot}/asterisk

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

if ( -d /usr/projects/jayenne/regress/asterisk/regression ) then
    echo "cd /usr/projects/jayenne/regress/asterisk/regression; svn update"
    cd /usr/projects/jayenne/regress/asterisk/regression
    svn update
else
    echo "mkdir -p /usr/projects/jayenne/regress/asterisk"
    mkdir -p /usr/projects/jayenne/regress/asterisk
    echo "svn co svn+ssh://ccscs8/ccs/codes/radtran/svn/asterisk/trunk/regression regression"
    svn co svn+ssh://ccscs8/ccs/codes/radtran/svn/asterisk/trunk/regression regression
endif

# Notes:
# ------------------------------------------------------------
# McKay uses the following:
# Purpose: Get Kerberos ticket
# Command: /usr/kerberos/bin/kinit -k -t /users/lmdm/.ssh/keytabfile lmdm@lanl.gov
# Output:  
