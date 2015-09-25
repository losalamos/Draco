# .bash_profile
#
# ~/.bash_profile is sourced only for login shells.
#------------------------------------------------------------------------------#

# Required values:
export DRACO_ENV_DIR=~/draco/environment

#------------------------------------------------------------------------------#
# CCS-2 standard setup
#
# - If this is an interactive shell then the environment variable $- should
#   contain an "i":
#------------------------------------------------------------------------------#
case ${-} in
  *i*)
    export INTERACTIVE=true
    ;;
  *) # Not an interactive shell
    export INTERACTIVE=false

    # Settings for non-interactive shells go here.  E.g.:
    # export PATH=/ccs/codes/radtran/vendors/subversion-1.8.13/bin:$PATH
    ;;
esac

#------------------------------------------------------------------------------#
# Load default Draco Environment
#------------------------------------------------------------------------------#
# module unload cmake ('module' may not be available for remote ssh connections)
if test "$INTERACTIVE" = true; then
    export prefered_term="konsole gnome-terminal xterm"
    if test -f ${DRACO_ENV_DIR}/bashrc/.bashrc; then
        source ${DRACO_ENV_DIR}/bashrc/.bashrc
    fi

    # Provide some bash functions
    if test -f ${DRACO_ENV_DIR}/bin/bash_functions.sh; then
        source ${DRACO_ENV_DIR}/bin/bash_functions.sh
    fi

    # User Customization goes here
    #
    # export EDITOR="emacs -nw"
    # export LPDEST=gumibears
    # alias xpdf='evince'

    # Set terminal title
    # echo -ne "\033]0;${nodename}\007"

    # prompt - see http://www.tldp.org/HOWTO/Bash-Prompt-HOWTO/
    #
    # if test "$TERM" = emacs || \
    #     test "$TERM" = dumb  || \
    #     test -z "`declare -f npwd | grep npwd`"; then
    #     export PS1="\h:\w [\!] % "
    #     export LS_COLORS=''
    # else
    #     found=`declare -f npwd | wc -l`
    #     if test ${found} != 0; then
    #         export PS1="\[\033[34m\]\h:\$(npwd) [\!] % \[\033[0m\]"
    #     fi
    # fi

    # Custom Modules
    #
    # module unuse /etc/modulefiles /usr/share/Modules/modulefiles
    # module swap cmake cmake/3.3.2
    # module list

fi
