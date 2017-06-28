#!/bin/bash
##-*- Mode: bash -*-
##---------------------------------------------------------------------------##
## File  : environment/bashrc/.bashrc
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##
##  Bash configuration file upon bash shell startup
##---------------------------------------------------------------------------##

#uncomment to debug this script.
#export verbose=true

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
   shopt -s cdspell # autocorrect spelling errors on cd command line.

   # Prevent creation of core files (ulimit -a to see all limits).
   ulimit -c 0

   ##------------------------------------------------------------------------##
   ## Common aliases
   ##------------------------------------------------------------------------##

   # Generic Settings

   alias ll='\ls -Flh'
   alias lt='\ls -Flth'
   alias ls='\ls -F'
   alias l.='\ls -h -d .*'

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
   alias xload="xload -label `hostname | sed -e 's/[.].*//'`"

   # Module related:
   alias ma='module avail'
   alias mls='module list'
   alias mld='module load'

   # Slurm related:
   # See http://hpcfaq.lanl.gov/faqs/view/58
   alias sqrun='squeue -o "%.7i %.8Q %.9P %.10u %.16a %.12j %.9B %.6C %.12L %S" -S "L" -t R'
   alias sqpend='squeue -o "%.7i %.8Q %.10u %.12a %.9P %.10j %.6C %.12S %.12L %.r" -t PD'
   alias pend='squeue -O "jobid:.7,prioritylong:.10,username:.10,account:.16,qos:.12,partition:.14,numnodes:.7,timelimit:.14,starttime:.22,reason:.20" -t PD'
   alias sqother='squeue -o "%.7i %.10u %.12T %.6C %10l %R" -S "-T" -t S,CG,CD,CF,CA,F,TO,PR,NF'
   alias sq='(echo RUNNING:;sqrun;echo "\nPENDING:";sqpend;echo "\nOther";sqother)'
   alias sqall='squeue -o "%.7i %.10u %.8T %.6C %.12L %S" -S "-t,-L"'
   alias showme='squeue -u ${USER} -o "%.7i %.10u %.8T %.6D %.6C %.12L %S" -S "-t,-L"'
   # checkjob analogs
   function chk() {
     job=$1
     if [ "$job" == "" ]; then
       echo No JobID Specified
     else
       scontrol show job $job
       #scontrol -d show job $job
       #scontrol -dd show job $job
     fi
   }
   # Show me all the possible qos and accounts
   function suser () {
     user=$1
     if [ "$user" == "" ]; then
       user=${USER}
     fi
     sacct -u $1 --format=JobID,Priority,JobName,Account,QOS,Partition,Start,Elapsed,State
   }
   function sshare_acct () {
     acct=$1
     if [ "$acct" == "" ]; then
       echo No Account Specified
     else
       sshare -m -a -A $1
     fi
   }
   function drmgroups () {
     user=$1
     if [ "$user" == "" ]; then
       user=${USER}
     fi
     sshare -m -u $user | grep $user | cut -f3 -d" " | uniq | sed '/^\s*$/d' | paste -d, -s
   }
   function sreason () {
     squeue -o "%.7i %.10u %.10T %r" -t PD,S,CG,CD,CF,CA,F,TO,PR,NF -j $1
   }
   function showuser () {
     user=$1
     if [ "$user" == "" ]; then
       user=${USER}
     fi
     sacctmgr list assoc user=$user format=Cluster,Account,Partition,QOS%-40
   }

   # Provide special ls commands if this is a color-xterm or compatible terminal.
   if test "${TERM}" != emacs &&
       test "${TERM}" != dumb; then
     # replace list aliases with ones that include colorized output.
     alias ll='\ls --color -Flh'
     alias l.='\ls --color -hd .*'
     alias lt='\ls --color -Flth'
     alias lt.='\ls --color -Flth .*'
     alias ls='\ls --color -F'
   fi

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
fi

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

  # Append PATHS (not linux specific, not ccs2 specific).
  add_to_path ${DRACO_ENV_DIR}/bin
  add_to_path ${DRACO_SRC_DIR}/tools

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
    tt-fey* | tt-login* | tr-fe* | tr-login* | nid0* )
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

  # Generic functions for loading/unloaded the default set of modules.
  fn_exists=`type dracoenv 2>/dev/null | head -n 1 | grep -c 'is a function'`
  if test $fn_exists = 0; then
    # only define if they do not already exist...
    function dracoenv ()
    {
      module load $dracomodules
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

##---------------------------------------------------------------------------##
## end of .bashrc
##---------------------------------------------------------------------------##
