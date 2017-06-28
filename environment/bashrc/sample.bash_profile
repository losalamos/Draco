# .bash_profile
#
# ~/.bash_profile is sourced only for login shells.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# CCS-2 standard setup
#
# - If this is an interactive shell then the environment variable $- should
#   contain an "i":
#------------------------------------------------------------------------------#
case ${-} in
  *i*) export INTERACTIVE=true ;;
  *)   export INTERACTIVE=false ;;
esac

#------------------------------------------------------------------------------#
# Draco Developer Environment
#------------------------------------------------------------------------------#
export DRACO_ENV_DIR=${HOME}/draco/environment
if [[ -f ${DRACO_ENV_DIR}/bashrc/.bashrc ]]; then
  source ${DRACO_ENV_DIR}/bashrc/.bashrc
fi

#------------------------------------------------------------------------------#
# User Customizations
#------------------------------------------------------------------------------#
# module unload cmake ('module' may not be available for remote ssh connections)
if test "$INTERACTIVE" = true; then

  # User Customization goes here

  # export EDITOR="emacs -nw"
  # export LPDEST=gumibears

  # Set terminal title
  # echo -ne "\033]0;${nodename}\007"

  # prompt - see http://www.tldp.org/HOWTO/Bash-Prompt-HOWTO/
  # if test "$TERM" = emacs || \
  #     test "$TERM" = dumb  || \
  #     test -z "`declare -f npwd | grep npwd`"; then
  #   export PS1="\h:\w [\!] % "
  #   export LS_COLORS=''
  # else
  #   found=`declare -f npwd | wc -l`
  #   if test ${found} != 0; then
  #     export PS1="\[\033[34m\]\h:\$(npwd) [\!] % \[\033[0m\]"
  #   fi
  # fi

  # Custom Modules
  # target="`uname -n | sed -e s/[.].*//`"
  # case ${target} in
  #   tt-fey* | tt-login*)
  #     module swap intel intel/15.0.3
  #     ;;
  #   *)
  #     module swap cmake cmake/3.4.0
  #     module list
  #     ;;
  # esac

  # LaTeX
  # export TEXINPUTS=$mydir:$TEXINPUTS
  # export BSTINPUTS=$mydir:$BSTINPUTS
  # export BIBINPUTS=$mydir:$BIBINPUTS

fi

#------------------------------------------------------------------------------#
# end ~/.bash_profile
#------------------------------------------------------------------------------#
