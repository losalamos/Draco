#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/push_repositories_xf.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# switch to group 'ccsrad' and set umask
if [[ $(id -gn) != ccsrad ]]; then
  exec sg ccsrad "$0 $*"
fi
umask 0007

# Locate the directory that this script is located in:
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Redirect all output to a log file.
timestamp=`date +%Y%m%d-%H%M`
target="`uname -n | sed -e s/[.].*//`"
logdir="$( cd $scriptdir/../../logs && pwd )"
logfile=$logdir/push_repositories_xf-$target-$timestamp.log
exec > $logfile
exec 2>&1

# import some bash functions
source $scriptdir/scripts/common.sh

echo -e "Executing $0 $*...\n"
echo "Group: `id -gn`"
echo -e "umask: `umask` \n"

#------------------------------------------------------------------------------#
# Q: How do I create a keytab that works with transfer 2.0
# A: See
#    https://rtt.lanl.gov/redmine/projects/draco/wiki/Kelly_Thompson#Generating-keytab-file-that-works-with-transfer-20
#    or see the comments at the end of this file.

# When kellyt runs this as a crontab, a special kerberos key must be used.
if [[ $USER == "kellyt" ]] ; then
    # Use a different cache location to avoid destroying any active user's
    # kerberos.
    export KRB5CCNAME=/tmp/regress_kerb_cache

    # Obtain kerberos authentication via keytab
    kinit -l 1h -kt $HOME/.ssh/xfkeytab transfer/${USER}push@lanl.gov
fi

# Sanity check
if test `klist -l | grep -c $USER` = 0; then
    die "You must have an active kerberos ticket to run this script."
fi

#------------------------------------------------------------------------------#
# Clone the github and gitlab repositories and push them to the red.

gitdir=/ccs/codes/radtran/git.${target}

if [[ -d $gitdir ]]; then
  run "cd $gitdir"
else
  die "GIT root directory not found. Expected to find $gitdir."
fi

repos="Draco.git jayenne.git capsaicin.git"

for repo in $repos; do

  # Remove the old tar file.
  if [[ -f $repo.tar ]]; then
    run "rm -f $repo.tar"
  fi

  if [[ -d $gitdir/$repo ]]; then

    # Tar it up
    run "tar -cvf $repo.tar $repo"

    # Ensure the new files have group rwX permissions.
    run "chgrp -R ccsrad $repo.tar"
    run "chmod -R g+rwX,o-rwX  $repo.tar"

    # Transfer the file via transfer.lanl.gov
    run "scp $repo.tar red@transfer.lanl.gov:"

  else
    echo "Warning: git mirror repository $gitdir/$repo was not found."
    echo "         Skipping to the next repository..."
  fi

done

echo -e "\n--------------------------------------------------------------------------------"
echo "All done."
echo "--------------------------------------------------------------------------------"


#------------------------------------------------------------------------------#
# Notes on using Transfer 2.0 (copied from the Draco wiki):
#
# * See http://transfer.lanl.gov.
# * See notes in cassio/Tools.rh/General/my_xf.pl

# <pre>
# % kinit -kt ~/.ssh/xfkeytab transfer/${USER}push@lanl.gov
# % ssh yellow@transfer.lanl.gov myfiles
# % scp /ccs/codes/radtran/svn/draco.kt.tar red@transfer.lanl.gov:
# </pre>

# * Using the transfer.sh script:

# <pre>
# draco/tools/transfer.sh push2r --recipients=kgt@lanl.gov -q `pwd`/<file>
# draco/tools/transfer.sh pullname <file>
# </pre>

# h3. Generating keytab file that works with transfer 2.0.

# * Part1: Obtain a business account via http://register.lanl.gov
# ** Under "Business Accounts Menu", click "Create a Business Account"
# ** Enter a name, (throughout this example, I will use <your_account>), then [Next].
# ** Enter Account Description, then [Next]
# *** [Confirm]
# ** Now, click "Manage Kerberos Principals"
# ** In the box, enter "transfer/" (without the quotes!) and [Create].
# ** At the "Are you sure..." prompt, [Create]
# ** You should have gotten a message that says,
# <pre>
#    Kerberos principal transfer/<your_account> was successfully Created.
# </pre>

# * Part 2: Register the principals
# ** Go to https://transfer.lanl.gov/principals
# ** Under "Add a New Principal", after the “transfer/”, enter your account name in the box:
# <pre>
#     transfer/<your_account>, click [Enroll Yellow Principal]
# </pre>
# ** If everything is correct, the four checkboxes below will get a green [Check]. If not, fix the first one that failed and try again. Repeat until success.

# * Part 3: Create the file
# ** Now, on the machine that you are going to use transfers*:
# ** *Note: Use ZNumber, not moniker here:*
# <pre>
#   % kadmin -p <ZNumber>@lanl.gov
#   Authenticating as principal <ZNumber>@lanl.gov with password.
#   Password for <ZNumber>@lanl.gov:  <cryptocard password>
# </pre>
# ** Now, at the kadmin: prompt, use the principal that you created, and append @lanl.gov (note in the xst command, “transfer_keytab” is the name of the file that will be created)
# <pre>
#    kadmin:   xst -k transfer_keytab transfer/<your_account>@lanl.gov
# Entry for principal transfer/<your_account>@lanl.gov with kvno 2, encryption type des-cbc-crc added to keytab WRFILE:transfer_keytab.
# Entry for principal transfer/<your_account>@lanl.gov with kvno 2, encryption type arcfour-hmac added to keytab WRFILE:transfer_keytab.
# Entry for principal transfer/<your_account>@lanl.gov with kvno 2, encryption type des3-cbc-sha1 added to keytab WRFILE:transfer_keytab.
# Entry for principal transfer/<your_account>@lanl.gov with kvno 2, encryption type aes128-cts-hmac-sha1-96 added to keytab
# WRFILE:transfer_keytab.
# Entry for principal transfer/<your_account>@lanl.gov with kvno 2, encryption type aes256-cts-hmac-sha1-96 added to keytab
# WRFILE:transfer_keytab.
#    kadmin: quit
# </pre>

# * Part 4: Test the file
# ** [first example, without the extracted key, no credentials]
# <pre>
#   % kdestroy
#   kdestroy: No credentials cache found while destroying cache

#   % ssh yellow@transfer.lanl.gov myfiles
#   yellow@transfer.lanl.gov's password:  <Ctl-C>
# </pre>
# ** [it did not work; now, use the keytab file]
# <pre>
#   % kinit -kt ./transfer_keytab transfer/<your_account>@lanl.gov

#   % ssh yellow@transfer.lanl.gov myfiles
#   ! No files found for zno 089990
# </pre>

# * Part 5: Request principal migration to Red
# ** https://transfer.lanl.gov/principal_request
# ** Now, log onto the red web page, https://transfer.lanl.gov/principals to register the same identity.

#------------------------------------------------------------------------------#
# End push_repositories_xf.sh
#------------------------------------------------------------------------------#
