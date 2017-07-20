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
   ;;
*) # Not an interactive shell
  export INTERACTIVE=false
  export DRACO_ENV_DIR=$HOME/draco/environment
  ;;
esac

#------------------------------------------------------------------------------#
# Draco developer environment
#------------------------------------------------------------------------------#
if test -f ${DRACO_ENV_DIR}/bashrc/.bashrc; then
  source ${DRACO_ENV_DIR}/bashrc/.bashrc
fi

#------------------------------------------------------------------------------#
# User customizations
#------------------------------------------------------------------------------#
if test "$INTERACTIVE" = true; then

  # Set terminal title
  # echo -ne "\033]0;${nodename}\007"

  # Aliases ---------------------------------------------------------------------#
  # alias xpdf='evince'
  # alias xload='xload -fg brown -fn 6x13 -geometry 180x100+1500+0'

  # drive your friends crazy...
  # alias vi='emacs -nw'

  # common cmake commands
  # alias cmakerel='cmake -DCMAKE_BUILD_TYPE=Release'
  # alias cmakerelfast='cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF'
  # alias cmakerwdi='cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo'
  # Use gcc's STL bounds checking
  # alias cmakebc='cmake -DGCC_ENABLE_GLIBCXX_DEBUG=ON'
  # Turn on extra debug info and floating point exception checking.
  # alias cmakefd='cmake -DDRACO_DIAGNOSTICS=7'

fi
