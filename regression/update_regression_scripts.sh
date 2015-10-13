#!/bin/bash

umask 0002

target="`uname -n | sed -e s/[.].*//`"
MYHOSTNAME="`uname -n`"
arch=`uname -m`

# Helper function
run () {
    echo $1
    if ! [ $dry_run ]; then eval $1; fi
}

# Ensure that the permissions are correct
case ${target} in
darwin-login*)
    echo "Please run regressions from darwin-fe instead of darwin-login."
    exit 1
    ;;
darwin-fe* | cn[0-9]*)
   /usr/projects/draco/vendors/keychain-2.7.1/keychain $HOME/.ssh/cmake_dsa
   if test -f $HOME/.keychain/$MYHOSTNAME-sh; then
      source $HOME/.keychain/$MYHOSTNAME-sh
   fi
   SVN=`which svn`
   # SVN=/projects/opt/subversion/1.8.10/bin/svn
   REGDIR=/usr/projects/draco/regress

   svnroot=/usr/projects/draco/regress/svn
   if ! test -d; then
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
      #echo "at svn+ssh://ccscs7/ccs/codes/radtran/svn/draco."
      #exit 1
      #fi
      #exit 0
      # cp pre-revprop-change.tmpl pre-revprop-change; chmod 775 \
      #    pre-revprop-change
      # vi pre-revprop-change --> comment out all code.
      # cd $svnroot
      # svnsync init file:///${svnroot}/jayenne svn+ssh://ccscs7/ccs/codes/radtran/svn/jayenne
      # svnsync sync file:///${svnroot}/jayenne
   fi

   run "${SVN}sync --non-interactive sync file://${svnroot}/draco"
   run "${SVN}sync --non-interactive sync file://${svnroot}/jayenne"
   run "${SVN}sync --non-interactive sync file://${svnroot}/capsaicin"
   # run "${SVN}sync --non-interactive sync file:///${svnroot}/asterisk"
   ;;
*)
   SVN=/ccs/codes/radtran/vendors/subversion-1.8.5/bin/svn
   REGDIR=/home/regress
   ;;
esac

# Update main regression scripts
run "cd ${REGDIR}/draco/config; ${SVN} update"
run "cd ${REGDIR}/draco/regression; ${SVN} update"
run "cd ${REGDIR}/draco/environment; ${SVN} update"
run "cd ${REGDIR}/jayenne/regression; ${SVN} update"
run "cd ${REGDIR}/capsaicin/scripts; ${SVN} update"
run "cd ${REGDIR}/asterisk/regression; ${SVN} update"
