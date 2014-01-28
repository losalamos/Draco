#!/bin/bash

# Is a kerberos ticket available?
# klist ....
#
#Ticket cache: FILE:/tmp/krb5cc_2017
#Default principal: kellyt@lanl.gov
#
#Valid starting     Expires            Service principal
#01/22/14 08:04:33  01/22/14 16:04:33  krbtgt/lanl.gov@lanl.gov
#	renew until 01/24/14 16:57:02

PUSH=/ccs/opt/x86_64/mercCmd-1.2.3/bin/push

# Helpful functions:
die () { echo "FATAL ERROR: $1"; exit 1;}

run () {
   echo $1
   if ! test $dry_run; then
      eval $1
   fi
}

# Working directory
if test -d /ccs/codes/radtran/svn; then
  run "cd /ccs/codes/radtran/svn"
else
  die "SVN root directory not found."
fi

if ! test -x $PUSH; then die "Cannot find Mercury executable push"; fi

# Repositories to push
repos="draco jayenne capsaicin"

#
for repo in $repos; do
   # Remove old hotcopy files/directories.
   if test -f ${repo}.hotcopy.tar; then
      run "rm -f ${repo}.hotcopy.tar"
   fi
   if test -d ${repo}.hotcopy; then
      run "rm -rf ${repo}.hotcopy"
   fi
   # Generate a repo hotcopy
   run "svnadmin hotcopy $repo ${repo}.hotcopy"
   # Tar it up and push via mercury.
   run "tar -cvf ${repo}.hotcopy.tar ${repo}.hotcopy"
   run "push id=${repo}.repo ${repo}.hotcopy.tar"
   # Ensure the new files have group rwX permissions.
   run "chgrp -R draco ${repo}.hotcopy.tar ${repo}.hotcopy"
   run "chmod -R g+rwX,o=g-w ${repo}.hotcopy.tar ${repo}.hotcopy"
done

