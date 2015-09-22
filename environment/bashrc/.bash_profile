###############################################################################
## .bash_profile
##
## $Id: .bash_profile 6776 2012-09-20 21:46:15Z kellyt $
###############################################################################

# This file is executed at login.
# .bashrc is executed at shell startup.

##---------------------------------------------------------------------------##
## INITIAL STUFF
##---------------------------------------------------------------------------##

# ~/.bash_profile: executed by bash(1) for login shells.

# If this is an interactive shell then the environment variable $- should
# contain an "i":
case ${-} in
*i*)
   export INTERACTIVE=true
   ;;
*) # Not an interactive shell
    export INTERACTIVE=false

   # Settings for non-interactive shells.  E.g.:
   # export PATH=/ccs/codes/radtran/vendors/subversion-1.8.13/bin:$PATH
   ;;
esac

#------------------------------------------------------------------------------#
# Load default Draco Environment
#------------------------------------------------------------------------------#
# module unload cmake ('module' may not be available for remote ssh connections)
if test "$INTERACTIVE" = true; then
    export prefered_term="konsole gnome-terminal xterm"
    export DRACO_ENV_DIR=~/draco/environment
    if test -f ${DRACO_ENV_DIR}/bashrc/.bashrc; then
        source ${DRACO_ENV_DIR}/bashrc/.bashrc
    fi

    # User Customization goes here
    #
    # export EDITOR="emacs -nw"
    #
    # alias xpdf='evince'
    #
    ##---------------------------------------------------------------------------##
    ## prompt - see http://www.tldp.org/HOWTO/Bash-Prompt-HOWTO/
    ##---------------------------------------------------------------------------##

    if test -f ${DRACO_ENV_DIR}/bin/bash_functions.sh; then
        source ${DRACO_ENV_DIR}/bin/bash_functions.sh
    fi

    if test "$TERM" = emacs || \
        test "$TERM" = dumb  || \
        test -z "`declare -f npwd | grep npwd`"; then
        export PS1="\h:\w [\!] % "
        export LS_COLORS=''
    else
        found=`declare -f npwd | wc -l`
        if test ${found} != 0; then
            export PS1="\[\033[34m\]\h:\$(npwd) [\!] % \[\033[0m\]"
        fi
    fi

fi
