#!/bin/bash

# Obtain kerberos authentication via keytab
kinit -f -l 1h -kt $HOME/.ssh/xfkeytab transfer/kellytpush@lanl.gov

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

# Repositories to push
repos="draco jayenne capsaicin"

#
for repo in $repos; do
   # Remove old hotcopy files/directories.
   # if test -f ${repo}.hotcopy.tar; then
   #    run "rm -f ${repo}.hotcopy.tar"
   # fi
   # if test -d ${repo}.hotcopy; then
   #    run "rm -rf ${repo}.hotcopy"
   # fi
   # # Generate a repo hotcopy
   # run "svnadmin hotcopy $repo ${repo}.hotcopy"
   # # Tar it up and push via mercury.
   # run "tar -cvf ${repo}.hotcopy.tar ${repo}.hotcopy"
   
   #run "${PUSH} push2r --quiet ${repo}.hotcopy.tar"
   run "scp ${repo}.hotcopy.tar red@transfer.lanl.gov:"

   # Ensure the new files have group rwX permissions.
   # run "chgrp -R draco ${repo}.hotcopy.tar ${repo}.hotcopy"
   # run "chmod -R g+rwX,o=g-w ${repo}.hotcopy.tar ${repo}.hotcopy"
done

