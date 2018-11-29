#!/bin/bash
##-*- Mode: bash -*-
##---------------------------------------------------------------------------##
## File  : environment/bashrc/.bashrc
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2018, Los Alamos National Security, LLC.
##         All rights are reserved.
##
## Bash configuration file upon bash shell startup
##
## Instructions (customization):
##
## 1. Setup
##    - Copy logic from draco/environment/bashrc/sample.bashrc and
##      draco/environment/bashrc/sample.bash_profile.
## 2. Override settings using the code found in the sample.bashrc.
##---------------------------------------------------------------------------##

#uncomment to debug this script.
#export verbose=true

##---------------------------------------------------------------------------##
## ENVIRONMENTS for interactive sessions
##---------------------------------------------------------------------------##

# If this is an interactive shell then the environment variable $- should
# contain an "i":
case ${-} in
  *i*)
    export INTERACTIVE=true
    if test -n "${verbose}"; then echo "in draco/environment/bashrc/.bashrc"; fi

    # Turn on checkwinsize
    shopt -s checkwinsize # autocorrect window size
    shopt -s cdspell      # autocorrect spelling errors on cd command line.
    shopt -s histappend   # append to the history file, don't overwrite it
    # shopt -s direxpand  # Doesn't work on toss22 machines. Move this command to
                          # .bashrc_toss3 and .bashrc_tt

    # don't put duplicate lines or lines starting with space in the history. See
    # bash(1) for more options
    HISTCONTROL=ignoreboth

    # for setting history length see HISTSIZE and HISTFILESIZE in bash(1)
    HISTSIZE=1000
    HISTFILESIZE=2000

    # Prevent creation of core files (ulimit -a to see all limits).
    # ulimit -c 0

    ##------------------------------------------------------------------------##
    ## Common aliases
    ##------------------------------------------------------------------------##

    # Generic Settings

    alias ll='\ls -Flh'
    alias lt='\ls -Flth'
    alias ls='\ls -F'
    alias la='\ls -A'
    alias l.='\ls -hd .*'
    alias lt.='ls -Flth .*'

    # alias a2ps='a2ps --sides=duplex --medium=letter'
    alias btar='tar --use-compress-program /usr/bin/bzip2'
    alias cpuinfo='cat /proc/cpuinfo'
    alias df='df -h'
    alias dirs='dirs -v'
    alias du='du -h --max-depth=1 --exclude=.snapshot'
    alias less='/usr/bin/less -r'
    alias mdstat='cat /proc/mdstat'
    alias meminfo='cat /proc/meminfo'
    alias mroe='more'
    nodename=`uname -n | sed -e 's/[.].*//g'`
    alias resettermtitle='echo -ne "\033]0;${nodename}\007"'

    # Module related:
    alias moduel='module'
    alias ma='module avail'
    alias mls='module list'
    alias mld='module load'
    alias mul='module unload'
    alias msh='module show'

    # set variable identifying the chroot you work in (used in the prompt below)
    if [ -z "${debian_chroot:-}" ] && [ -r /etc/debian_chroot ]; then
      debian_chroot=$(cat /etc/debian_chroot)
    fi

    # If this is an xterm set the title to user@host:dir
    case "$TERM" in
      xterm*|rxvt*)
        # PS1="\[\e]0;${debian_chroot:+($debian_chroot)}\u@\h: \w\a\]$PS1"
        echo -ne "\033]0;${nodename}\007"
        ;;
      *)
        ;;
    esac

    # Provide special ls commands if this is a color-xterm or compatible
    # terminal.

    # 1. Does the current terminal support color?
    if [ -x /usr/bin/tput ] && tput setaf 1 >&/dev/null; then
      # We have color support; assume it's compliant with Ecma-48
      # (ISO/IEC-6429). (Lack of such support is extremely rare, and such a case
      # would tend to support setf rather than setaf.)
      color_prompt=yes
    fi

    # 2. Override color_prompt for special values of $TERM
    case "$TERM" in
      xterm-color|*-256color) color_prompt=yes;;
      emacs|dumb)
        color_prompt=no
        LS_COLORS=''
        ;;
    esac

    # if ! [ -x /usr/bin/dircolors ]; then
    #   color_prompt=no
    # fi

    if [[ "${color_prompt:-no}" == "yes" ]]; then

      # Use custom colors if provided.
      test -r ~/.dircolors && eval "$(dircolors -b ~/.dircolors)" || eval "$(dircolors -b)"

      # append --color option to some aliased commands

      alias ll='\ls -Flh --color'
      alias lt='\ls -Flth --color'
      alias ls='\ls -F --color'
      alias la='\ls -A --color'
      alias l.='\ls -hd --color .*'
      alias lt.='ls -Flth --color .*'

      alias grep='grep --color=auto'
      alias fgrep='fgrep --color=auto'
      alias egrep='egrep --color=auto'

      # colored GCC warnings and errors
      export GCC_COLORS='error=01;31:warning=01;35:note=01;36:caret=01;32:locus=01:quote=01'

      # Colorized prompt (might need some extra debian_chroot stuff -- see wls
      # example).
      if [ -z "${debian_chroot:-}" ] && [ -r /etc/debian_chroot ]; then
        debian_chroot=$(cat /etc/debian_chroot)
      fi

      if [ "$color_prompt" = yes ]; then
        PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ '
      else
        PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w\$ '
      fi

    fi
    unset color_prompt

    ;; # end case 'interactive'

  ##---------------------------------------------------------------------------##
  ## ENVIRONMENTS for non interactive sessions
  ##---------------------------------------------------------------------------##

  *) # Not an interactive shell (e.g. A PBS shell?)
    export INTERACTIVE=false
    ;;
esac

##---------------------------------------------------------------------------##
## ENVIRONMENTS - bash functions (all interactive sessions)
##---------------------------------------------------------------------------##

# Bash functions are not inherited by subshells.
if [[ ${INTERACTIVE} ]]; then

  # Attempt to find DRACO
  if ! [[ $DRACO_SRC_DIR ]]; then
    _BINDIR=`dirname "$BASH_ARGV"`
    export DRACO_SRC_DIR=`(cd $_BINDIR/../..;pwd)`
    export DRACO_ENV_DIR=${DRACO_SRC_DIR}/environment
  fi

  # Common bash functions and alias definitions
  source ${DRACO_ENV_DIR}/bin/bash_functions.sh
  source ${DRACO_ENV_DIR}/../regression/scripts/common.sh

  # aliases and bash functions for working with slurm
  if !  [[ `which squeue 2>&1 | grep -c "no squeue"` == 1 ]] &&
    [[ `which squeue | grep -c squeue` -gt 0 ]]; then
    source ${DRACO_ENV_DIR}/bashrc/.bashrc_slurm
  fi
fi

##---------------------------------------------------------------------------##
## ENVIRONMENTS - once per login
##---------------------------------------------------------------------------##

# Darwin salloc inherits the user environment, so we need to bypass the
# "already-done" logic
if [[ ${SLURM_CLUSTER_NAME} == "darwin" ]]; then
  export DRACO_BASHRC_DONE=no
fi

if [[ ${DRACO_BASHRC_DONE:-no} == no ]] && [[ ${INTERACTIVE} == true ]]; then

  # Clean up the default path to remove duplicates
  tmpifs=$IFS
  oldpath=$PATH
  export PATH=/bin
  IFS=:
  for dir in $oldpath; do
    if test -z "`echo $PATH | grep $dir`" && test -d $dir; then
      export PATH=$PATH:$dir
    fi
  done
  IFS=$tmpifs
  unset tmpifs
  unset oldpath

  # Avoid double colon in PATH
  export PATH=`echo ${PATH} | sed -e 's/[:]$//'`
  export LD_LIBRARY_PATH=`echo ${LD_LIBRARY_PATH} | sed -e 's/[:]$//'`

  # Append PATHS (not linux specific, not ccs2 specific).
  add_to_path ${DRACO_ENV_DIR}/bin
  add_to_path ${DRACO_SRC_DIR}/tools

  # Tell wget to use LANL's www proxy (see
  # trac.lanl.gov/cgi-bin/ctn/trac.cgi/wiki/SelfHelpCenter/ProxyUsage)
  # export http_proxy=http://wpad.lanl.gov/wpad.dat
  current_domain=`awk '/^domain/ {print $2}' /etc/resolv.conf`
  #  found=`nslookup proxyout.lanl.gov | grep -c Name`
  #  if test ${found} == 1; then
  if [[ ${current_domain} == "lanl.gov" ]]; then
    export http_proxy=http://proxyout.lanl.gov:8080
    export https_proxy=$http_proxy
    export HTTP_PROXY=$http_proxy
    export HTTPS_PROXY=$http_proxy
    export no_proxy=".lanl.gov"
    export NO_PROXY=$no_proxy
  fi

  # cd paths - disable here, let user choose in ~/.bashrc
  CDPATH=

  # Hooks for clang-format as git commit hook:
  # Possible values: ON, TRUE, OFF, FALSE, DIFF (the default value is ON).
  export DRACO_AUTO_CLANG_FORMAT=ON

  # Silence warnings from GTK/Gnome
  export NO_AT_BRIDGE=1

  ##---------------------------------------------------------------------------##
  ## ENVIRONMENTS - machine specific settings
  ##---------------------------------------------------------------------------##
  target="`uname -n | sed -e s/[.].*//`"
  arch=`uname -m`

  case ${target} in

    # machine with GPUs
    # backend nodes with GPUs are cn[1-4].
    darwin-fe* | cn[0-9]*)
      if [[ $verbose ]]; then echo "this is Darwin"; fi
      source ${DRACO_ENV_DIR}/bashrc/.bashrc_darwin_fe
      ;;

    # Pinto | Wolf
    pi* | wf* | lu* )
      source ${DRACO_ENV_DIR}/bashrc/.bashrc_toss22
      ;;

    # Badger | Fire | Grizzly | Ice | Snow
    ba* | cy* | fi* | gr* | ic* | sn* )
      source ${DRACO_ENV_DIR}/bashrc/.bashrc_toss3
      ;;

    # wtrw and rfta
    red-wtrw* | rfta* | redcap* )
      source ${DRACO_ENV_DIR}/bashrc/.bashrc_rfta
      ;;
    # trinitite (tt-fey) | trinity (tr-fe)
    tt-fey* | tt-login* | tr-fe* | tr-login* | nid* )
      source ${DRACO_ENV_DIR}/bashrc/.bashrc_tt
      ;;
    # rzuseq
    rzuseq*)
      source ${DRACO_ENV_DIR}/bashrc/.bashrc_bgq
      ;;

    # Assume CCS machine (ccscs[0-9] or personal workstation)
    *)
      if [[ -d /ccs/codes/radtran ]]; then
        # assume this is a CCS LAN machine (64-bit)
        if test `uname -m` = 'x86_64'; then
          # draco environment only supports 64-bit linux...
          source ${DRACO_ENV_DIR}/bashrc/.bashrc_linux64
        else
          echo "Draco's environment is not fully supported on 32-bit Linux."
          echo "Module support may not be available. Email kgt@lanl.gov for more information."
          # source ${DRACO_ENV_DIR}/bashrc/.bashrc_linux32
        fi
      elif [[ -d /usr/projects/draco ]]; then
        # XCP machine like 'toolbox'?
        source ${DRACO_ENV_DIR}/bashrc/.bashrc_linux64
      fi
      export NoModules=1
      ;;

  esac

  source ${DRACO_ENV_DIR}/bashrc/bash_functions2.sh
  dracoenv

  # Mark that we have already done this setup
  export DRACO_BASHRC_DONE=yes

fi

# provide some bash functions (dracoenv, rmdracoenv) for non-interactive
# sessions.
source ${DRACO_ENV_DIR}/bashrc/bash_functions2.sh


if test -n "${verbose}"; then echo "done with draco/environment/bashrc/.bashrc"; fi

##---------------------------------------------------------------------------##
## end of .bashrc
##---------------------------------------------------------------------------##
