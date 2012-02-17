#!/bin/tcsh -f

# See jayenne-project/regression/sync_repository.csh.

# setenv MYHOSTNAME `uname -n`

# # Credentials via Keychain (SSH)
# /usr/projects/draco/vendors/keychain-2.7.1/keychain $HOME/.ssh/cmake_dsa
# source $HOME/.keychain/$MYHOSTNAME-csh

# # Setup directory structure
# if ( -d /usr/projects/jayenne/regress/cvsroot ) then
#    :
# else
#    echo "install -d /usr/projects/jayenne/regress/cvsroot"
#    install -d /usr/projects/jayenne/regress/cvsroot
# endif

# cd  /usr/projects/jayenne/regress/cvsroot
# # scp -r ccscs8:/ccs/codes/radtran/cvsroot/CVSROOT .
# # scp -r ccscs8:/ccs/codes/radtran/cvsroot/LOG .
# # scp -r ccscs8:/ccs/codes/radtran/cvsroot/bin .
# # scp -r ccscs8:/ccs/codes/radtran/cvsroot/draco .
# # scp -r ccscs8:/ccs/codes/radtran/cvsroot/clubimc .
# # scp -r ccscs8:/ccs/codes/radtran/cvsroot/wedgehog .
# # scp -r ccscs8:/ccs/codes/radtran/cvsroot/milagro .

# setenv packages "CVSROOT LOG bin draco clubimc wedgehog milagro imcdoc"
# foreach pkg  ( $packages )
#    echo "----------------------------------------------------------------------"
#    echo "rsync -avzr ccscs8:/ccs/codes/radtran/cvsroot/$pkg/ $pkg"
#    rsync -avzr --delete ccscs8:/ccs/codes/radtran/cvsroot/$pkg/ $pkg
# end


# # McKay uses the following:
# # Purpose: Get Kerberos ticket
# # Command: /usr/kerberos/bin/kinit -k -t /users/lmdm/.ssh/keytabfile lmdm@lanl.gov
# # Output:  
