# .bashrc
#
# ~/.bashrc is sourced only for non-login shells.
#-------------------------------------------------------------------------------
if test -f /etc/bashrc; then
  source /etc/bashrc
fi

#------------------------------------------------------------------------------#
# CCS-2 standard setup
#
# - If this is an interactive shell then the environment variable $- should
#   contain an "i":
#------------------------------------------------------------------------------#
case ${-} in
*i*)
   export INTERACTIVE=true
   # If this is an interactive shell and DRACO_ENV_DIR isn't set. Assume that we
   # need to source the .bash_profile.
   if test -z "${DRACO_ENV_DIR}" && test -f ${HOME}/.bash_profile; then
       source $HOME/.bash_profile
   fi
   # Ensure Draco's bash_functions are declared in the non-login subshell.
   if test -f ${DRACO_ENV_DIR}/bin/bash_functions.sh; then
       source ${DRACO_ENV_DIR}/bin/bash_functions.sh
   fi
   ;;
*) # Not an interactive shell
    export INTERACTIVE=false
   ;;
esac

#------------------------------------------------------------------------------#
# User customizations
#------------------------------------------------------------------------------#
if test "$INTERACTIVE" = true; then

    # Set terminal title
    # echo -ne "\033]0;${nodename}\007"

fi
