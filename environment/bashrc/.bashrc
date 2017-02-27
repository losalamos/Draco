#!/bin/bash
##-*- Mode: bash -*-
##---------------------------------------------------------------------------##
## File  : environment/bashrc/.bashrc
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016, Los Alamos National Security, LLC.
##         All rights are reserved.
##
##  Bash configuration file upon bash shell startup
##---------------------------------------------------------------------------##

## Instructions (customization)
##
## Before sourcing this file, you may wish to set the following
## variables to customize your environment (ie: set in ~/.bashrc
## before sourcing this file).
##
## $prefered_term - a space delimited list of terminal names.  The
##           default list is "gnome-terminal knosole xterm".  You can
##           modify the order of this list or remove items.  If you
##           add a new terminal you will need to modify this file to
##           set the optional parameters.

##---------------------------------------------------------------------------##
## ENVIRONMENTS for interactive sessions
##---------------------------------------------------------------------------##

# If this is an interactive shell then the environment variable $-
# should contain an "i":
case ${-} in
*i*)
   export INTERACTIVE=true
   if test -n "${verbose}"; then
      echo "in draco/environment/bashrc/.bashrc"
   fi

   # Turn on checkwinsize
   shopt -s checkwinsize # autocorrect window size
   #shopt -s cdspell # autocorrect spelling errors on cd command line.

   # Prevent creation of core files (ulimit -a to see all limits).
   ulimit -c 0

   ##------------------------------------------------------------------------##
   ## Common aliases
   ##------------------------------------------------------------------------##

   # Generic Settings

   alias ll='ls -Flh'
   alias lt='ls -Flth'
   alias ls='ls -F'
   alias l.='ls -h -d .*'

   # alias a2ps='a2ps --sides=duplex --medium=letter'
   alias btar='tar --use-compress-program /usr/bin/bzip2'
   alias cd..='cd ..'
   alias cpuinfo='cat /proc/cpuinfo'
   alias df='df -h'
   alias dirs='dirs -v'
   alias dmesg='dmesg -s 65536'
   alias du='du -h --max-depth=1 --exclude=.snapshot'
   alias free='free -m'
   alias hosts='cat /etc/hosts'
   alias hpss='echo Try using psi instead of hpss'
   alias ldmesg='dmesg -s 65536 | less'
   alias less='/usr/bin/less -r'
   alias mdstat='cat /proc/mdstat'
   alias meminfo='cat /proc/meminfo'
   alias mroe='more'
   alias print='lp'
   nodename=`uname -n | sed -e 's/[.].*//g'`
   alias resettermtitle='echo -ne "\033]0;${nodename}\007"'
   alias sdl='DISPLAY=127.0.0.1:0.0;echo The DISPLAY value is now: $DISPLAY'
   alias xload="xload -label `hostname | sed -e 's/[.].*//'`"

   # Provide special ls commands if this is a color-xterm or compatible terminal.
   if test "${TERM}" != emacs &&
       test "${TERM}" != dumb; then
     # replace list aliases with ones that include colorized output.
     alias ll='ls --color -Flh'
     alias l.='ls --color -hd .*'
     alias lt='ls --color -Flth'
     alias lt.='ls --color -Flth .*'
     alias ls='ls --color -F'
   fi

   # Turquoise network
   #alias mapache='ssh -t -X wtrw.lanl.gov ssh mp-fe1'
   #alias tscp='scp $1 turq-fta1.lanl.gov:/scratch/$USERNAME/$1'
   #alias trsync='rsync -avz -e ssh --protocol=20 $1 turq-fta1.lanl.gov:/scratch/$USERNAME/$1'

   ;; # end case 'interactive'

##---------------------------------------------------------------------------##
## ENVIRONMENTS for non interactive sessions
##---------------------------------------------------------------------------##

*) # Not an interactive shell (e.g. A PBS shell?)
   export INTERACTIVE=false
   # return
   ;;
esac
#export verbose=true

##---------------------------------------------------------------------------##
## ENVIRONMENTS - once per login
##---------------------------------------------------------------------------##

if test ${DRACO_BASHRC_DONE:-no} = no && test ${INTERACTIVE} = true; then

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

  # Attempt to find DRACO
  if test -z "$DRACO_SRC_DIR"; then
    _BINDIR=`dirname "$BASH_ARGV"`
    export DRACO_SRC_DIR=`(cd $_BINDIR/../..;pwd)`
  fi
  if test -z "$DRACO_ENV_DIR"; then
    export DRACO_ENV_DIR=${DRACO_SRC_DIR}/environment
  fi

  # Append PATHS (not linux specific, not ccs2 specific).
  extradirs="${DRACO_ENV_DIR}/bin ${DRACO_SRC_DIR}/tools"
  for mydir in ${extradirs}; do
    if test -z "`echo $PATH | grep $mydir`" && test -d $mydir; then
      export PATH=${PATH}:${mydir}
    fi
  done

  # Tell wget to use LANL's www proxy (see
  # trac.lanl.gov/cgi-bin/ctn/trac.cgi/wiki/SelfHelpCenter/ProxyUsage)
  # export http_proxy=http://wpad.lanl.gov/wpad.dat
  found=`nslookup proxyout.lanl.gov | grep -c Name`
  if test ${found} == 1; then
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

  ##---------------------------------------------------------------------------##
  ## ENVIRONMENTS - machine specific settings
  ##---------------------------------------------------------------------------##
  target="`uname -n | sed -e s/[.].*//`"
  arch=`uname -m`

  case ${target} in

    # machine with GPUs
    # backend nodes with GPUs are cn[1-4].
    darwin-fe* | cn[0-9]*)
      source ${DRACO_ENV_DIR}/bashrc/.bashrc_darwin_fe
      ;;

    # Mapache | Moonlight | Mustang | Pinto | Wolf | Luna
    lu* | ml* | pi* | wf* )
      source ${DRACO_ENV_DIR}/bashrc/.bashrc_toss22
      ;;

    # Snow | Fire | Ice
    sn* | fi* | ic* )
      source ${DRACO_ENV_DIR}/bashrc/.bashrc_toss3
      ;;

    # wtrw and rfta
    red-wtrw* | rfta* | redcap* )
      source ${DRACO_ENV_DIR}/bashrc/.bashrc_rfta
      ;;
    # trinitite (tt-fey) | trinity (tr-fe)
    tt-fey* | tt-login* | tr-fe* | tr-login*)
      source ${DRACO_ENV_DIR}/bashrc/.bashrc_tt
      ;;
    # rzuseq
    rzuseq*)
      source ${DRACO_ENV_DIR}/bashrc/.bashrc_bgq
      ;;

    # Assume CCS machine (ccscs[0-9] or personal workstation)
    *)
      if test -d /ccs/codes/radtran; then
        # assume this is a CCS LAN machine (64-bit)
        if test `uname -m` = 'x86_64'; then
          # draco environment only supports 64-bit linux...
          source ${DRACO_ENV_DIR}/bashrc/.bashrc_linux64
        else
          echo "Draco's environment is not fully supported on 32-bit Linux."
          echo "Module support may not be available. Email kgt@lanl.gov for more information."
          source ${DRACO_ENV_DIR}/bashrc/.bashrc_linux32
        fi
      fi
      export NoModules=1
      ;;

  esac

  # Generic functions for loading/unloaded the default set of modules.
  fn_exists=`type dracoenv 2>/dev/null | head -n 1 | grep -c 'is a function'`
  if test $fn_exists = 0; then
    # only define if they do not already exist...
    function dracoenv ()
    {
      for m in $dracomodules; do
        module load $m
      done
    }
    function rmdracoenv ()
    {
      # unload in reverse order.
      mods=( ${dracomodules} )
      for ((i=${#mods[@]}-1; i>=0; i--)); do
        loaded=`echo $LOADEDMODULES | grep -c ${mods[$i]}`
        if test $loaded = 1; then
          module unload ${mods[$i]}
        fi
      done
    }
  fi
  dracoenv

  # Mark that we have already done this setup
  export DRACO_BASHRC_DONE=yes

fi

# Common bash functions and alias definitions
source ${DRACO_ENV_DIR}/bin/bash_functions.sh

##---------------------------------------------------------------------------##
## end of .bashrc
##---------------------------------------------------------------------------##
