# .bashrc_functions
#
# .bashrc_functions is sourced by interactive shells from
# .bash_profile and .bashrc
#------------------------------------------------------------------------------#

# My Functions ----------------------------------------------------------------#

# function eapde ()
# {
#   export EAP_INHIBIT_KSM=1
#   export L_EAP_ENV=new
#   source /usr/projects/eap/dotfiles/.bashrc
#   module unload paraview
#   module list
# }
# export -f eapde

#------------------------------------------------------------------------------#
# function dos2unix ()
# {
#   dosfile=$1
#   unixfile=$dosfile
#   tmp=`mktemp`
#   awk '{ sub("\r$", ""); print }' $dosfile > $tmp
#   cp $tmp $unixfile
# }
# export -f dos2unxi

# function unix2dos ()
# {
#   unixfile=$1
#   dosfile=$unixfile
#   tmp=`mktemp`
#   awk 'sub("$", "\r")' $unixfile > $tmp
#   cp $tmp $dosfile
# }
# export -f unix2dos

#------------------------------------------------------------------------------#
# function keychain()
# {
#   # https://www.funtoo.org/Keychain
#   if [[ -d $VENDOR_DIR/keychain-2.8.5 ]]; then
#     eval "$VENDOR_DIR/keychain-2.8.5/keychain --eval --agents ssh $HOME/.ssh/id_rsa"
#   fi
# }
# export -f keychain

#------------------------------------------------------------------------------#
# End .bash_functions
#------------------------------------------------------------------------------#
