###############################################################################
## .bash_profile
##
## $Id$
###############################################################################

# This file is executed at login.
# .bashrc is executed at shell startup.

##---------------------------------------------------------------------------##
## INITIAL STUFF
##---------------------------------------------------------------------------##

# ~/.bash_profile: executed by bash(1) for login shells.

# If this is an interactive shell then the environment variable $-
# should contain an "i":
case ${-} in 
*i*)
   export INTERACTIVE=true
   verbose=
   if test -n "${verbose}"; then
      echo "in .bash_profile"
   fi
   ;;
*) # Not an interactive shell
   export INTERACTIVE=false
   ;;
esac

##---------------------------------------------------------------------------##
## Establish environment for login shell
##---------------------------------------------------------------------------##
if [ -f ~/draco/environment/bashrc/.bashrc ]; then
  source ~/draco/environment/bashrc/.bashrc
fi
