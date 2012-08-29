##-*- Mode: bash -*-
##---------------------------------------------------------------------------##
## .bashrc - my bash configuration file upon bash shell startup
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
## ENVIRONMENTS for all sessions
##---------------------------------------------------------------------------##

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

# Append PATHS (not linux specific, not ccs2 specific).
extradirs="${DRACO_SRC_DIR}/environment/bin ${DRACO_SRC_DIR}/tools /usr/X11R6/bin /usr/lanl/bin"
for mydir in ${extradirs}; do
   if test -z "`echo $PATH | grep $mydir`" && test -d $mydir; then
      export PATH=${PATH}:${mydir}
   fi
done

# Tell wget to use LANL's www proxy (see
# trac.lanl.gov/cgi-bin/ctn/trac.cgi/wiki/SelfHelpCenter/ProxyUsage) 
# export http_proxy=http://wpad.lanl.gov/wpad.dat
export http_proxy=http://proxyout.lanl.gov:8080
export https_proxy=$http_proxy
export HTTP_PROXY=$http_proxy
export HTTPS_PROXY=$http_proxy
export http_no_proxy="*.lanl.gov"
export no_proxy=lanl.gov
export NO_PROXY=$no_proxy
# See help page for how to setup subversion

# cd paths - disable here, let user choose in ~/.bashrc
CDPATH=

##---------------------------------------------------------------------------##
## ENVIRONMENTS for interactive sessions
##---------------------------------------------------------------------------##

# If this is an interactive shell then the environment variable $-
# should contain an "i":
case ${-} in 
*i*)
   export INTERACTIVE=true
   export verbose=
   if test -n "${verbose}"; then
      echo "in .bashrc"
   fi

   # Turn on checkwinsize
   shopt -s checkwinsize # autocorrect window size
   #shopt -s cdspell # autocorrect spelling errors on cd command line.

   # X server resources
   if test -f ${HOME}/.Xdefaults; then
       if test -x /usr/X11R6/bin/xrdb; then
           if test ! "${DISPLAY}x" = "x"; then
               /usr/X11R6/bin/xrdb ${HOME}/.Xdefaults
           fi
       fi
   fi

   ##---------------------------------------------------------------------------##
   ## aliases
   ##---------------------------------------------------------------------------##

   # Generic Settings

   alias ll='ls -Fl'
   alias lt='ls -Flt'
   alias ls='ls -F'
   alias l.='ls -h -d .*'

   alias a2ps='a2ps --sides=duplex --medium=letter'
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
   alias watchioblocks='ps -eo stat,pid,user,command | egrep "^STAT|^D|^R"'
   #alias whatsmyip='wget -O - -q myip.dk | grep Box | grep div | egrep -o [0-9.]+'
   alias wmdstat='watch -n 2 "cat /proc/mdstat"'
   alias xload="xload -label `hostname | sed -e 's/[.].*//'`"

   if test -x /ccs/codes/marmot/magicdraw/MagicDraw_UML_9.0/bin/mduml; then
       alias magicdraw='/ccs/codes/marmot/magicdraw/MagicDraw_UML_9.0/bin/mduml'
   fi
   if test -x /usr/bin/kghostview; then
       alias gv='/usr/bin/kghostview'
       alias ghostview='/usr/bin/kghostview'
   fi
   if test -x /opt/bin/gstat; then
       alias wgstat='watch /opt/bin/gstat -1a'
       alias linkstat='gsh dmesg | grep "Link is up at" | sort -u'
       alias interactive_qsub='xterm -e qsub -I -l nodes=1:ppn=4 &'
       alias cqs='echo -e "\nCurrent Queuing System: $PREFERED_QUEUE_SYSTEM \n"'
   fi

   # Provide special ls commands if this is a color-xterm or compatible terminal.
   if test "${TERM}" != emacs && 
       test "${TERM}" != dumb; then
   # replace list aliases with ones that include colorized output.
       alias ll='ls --color -Fl'
       alias l.='ls --color -aFl'
       alias lt='ls --color -Flt'
       alias lt.='ls --color -aFlt'
       alias ls='ls --color -F'
   fi

   source ${DRACO_SRC_DIR}/environment/bin/bash_functions.sh

   # Aliases for machines

   # No need to use ssh to pop a terminal from the current machine
   # alias ${target}='${term} ${term_opts}'

   # Turquise network
   alias mapache='ssh -t -X wtrw.lanl.gov ssh mp-fe1'
   alias tscp='scp $1 turq-fta1.lanl.gov:/scratch/$USERNAME/$1'
   alias trsync='rsync -avz -e ssh --protocol=20 $1 turq-fta1.lanl.gov:/scratch/$USERNAME/$1'

##---------------------------------------------------------------------------##
## prompt - see http://www.tldp.org/HOWTO/Bash-Prompt-HOWTO/
##---------------------------------------------------------------------------##

   if test "$TERM" = emacs || \
       test "$TERM" = dumb  || \
       test -z "`declare -f npwd | grep npwd`"; then
       export PS1="\h:\w [\!] % "
       export LS_COLORS=''
   else
       export PS1="\[\033[34m\]\h:\$(npwd) [\!] % \[\033[0m\]"
   fi
   ;;

##---------------------------------------------------------------------------##
## ENVIRONMENTS for non interactive sessions
##---------------------------------------------------------------------------##

*) # Not an interactive shell (e.g. A PBS shell?)
   export INTERACTIVE=false
   # return
   ;;
esac

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
## Machine specific settings
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

target="`uname -n | sed -e s/[.].*//`"
arch=`uname -m`

case ${target} in

# machine with GPUs
# backend nodes with GPUs are cn[1-4].
darwin | cn[0-9]*)
   source ${DRACO_SRC_DIR}/environment/bashrc/.bashrc_darwin
   ;; 

# RoadRunner machines
rra[0-9]*a)
    source ${DRACO_SRC_DIR}/environment/bashrc/.bashrc_rr
    ;;
rr-dev-fe)
    source ${DRACO_SRC_DIR}/environment/bashrc/.bashrc_rr_dev
    ;;

#Yellowrail machines
yr-fe1 | yra[0-9]* )
    source ${DRACO_SRC_DIR}/environment/bashrc/.bashrc_yr
    ;;

#TLCC machines
tu-fe* | tua[0-9]* | hu-fe[1-2] | hu*[0-9]* | ty-fe* | ty[0-9]*)
    source ${DRACO_SRC_DIR}/environment/bashrc/.bashrc_tlcc
    ;;

# Cielito
ct-fe[0-9] | ct-login[0-9])
    source ${DRACO_SRC_DIR}/environment/bashrc/.bashrc_ct
    ;;

# Luna
lu-fe[0-9] | lua[0-9]* | ml-fey | ml*)
    source ${DRACO_SRC_DIR}/environment/bashrc/.bashrc_lu
    ;;

# Mapache

# Assume CCS machine (ccscs[0-9] or personal workstation)
*)
    if test -d /ccs/codes/radtran; then 
        # assume this is a CCS LAN machine (64-bit)
        if test `uname -m` = 'x86_64'; then
          # draco environment only supports 64-bit linux...
          source ${DRACO_SRC_DIR}/environment/bashrc/.bashrc_linux64
        else
          echo "Draco's environment no longer supports 32-bit Linux."
          echo "Module support not available. Email kgt@lanl.gov for more information."
        fi
    fi
    ;;

esac

# Only print the loaded modules if this is an interactive session.
case ${-} in 
*i*)
    if test -n "$MODULESHOME"; then
      module list
    fi
    ;;
esac

##---------------------------------------------------------------------------##
## end of .bashrc
##---------------------------------------------------------------------------##
