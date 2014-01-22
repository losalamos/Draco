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

# Working directory
if test -d /ccs/codes/radtran/svn; then
  cd /ccs/codes/radtran/svn
else
  echo "FATAL ERROR.  SVN root directory not found."
  exit 1
fi

# Repositories to push
repos="draco jayenne capsaicin"

#
for repo in $repos; do
   # Remove old hotcopy files/directories.
   if test -f ${repo}.hotcopy.tar; then
      rm -f ${repo}.hotcopy.tar
   fi
   if test -d ${repo}.hotcopy; then
      rm -rf ${repo}.hotcopy
   fi
   # Generate a repo hotcopy
   svnadmin hotcopy $repo ${repo}.hotcopy
   # Tar it up and push via mercury.
   tar -cvf ${repo}.hotcopy.tar ${repo}.hotcopy
   push id=${repo}.repo ${repo}.hotcopy.tar
   # Ensure the new files have group rwX permissions.
   chgrp -R ${repo}.hotcopy.tar ${repo}.hotcopy
   chmod -R g+rwX,o=g-w ${repo}.hotcopy.tar ${repo}.hotcopy
done

