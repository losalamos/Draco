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
## If this is a non interactive shell then don't process this rc file
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
   ;;
*) # Not an interactive shell (e.g. A PBS shell?)
   export INTERACTIVE=false
   return
   ;;
esac

##---------------------------------------------------------------------------##
## ENVIRONMENTS
##---------------------------------------------------------------------------##

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

# Avoid double colon in PATH
export PATH=`echo ${PATH} | sed -e 's/[:]$//'`
export LD_LIBRARY_PATH=`echo ${LD_LIBRARY_PATH} | sed -e 's/[:]$//'`

# Append PATHS (not linux specific, not ccs2 specific).

if test -z "$DRACO_SRC_DIR"; then
  _BINDIR=`dirname "$BASH_ARGV"`
  export DRACO_SRC_DIR=`(cd $_BINDIR/../..;pwd)`
fi

extradirs="${DRACO_SRC_DIR}/environment/bin /usr/X11R6/bin /usr/lanl/bin"
for mydir in ${extradirs}; do
   if test -z "`echo $PATH | grep $mydir`" && test -d $mydir; then
      export PATH=${PATH}:${mydir}
   fi
done

# set variable with my moniker.
export USERNAME=`basename ${HOME}`

# Remove all permissions for world and group for files I create.
umask 077

# Tell the Draco build system to use all available cores when
# compiling.
if test -f /proc/cpuinfo; then
  export nj=`cat /proc/cpuinfo | grep processor | wc -l`
fi

# Tell wget to use LANL's www proxy
# export http_proxy=http://wpad.lanl.gov/wpad.dat
# export http_proxy=http://proxyout.lanl.gov:8080

##---------------------------------------------------------------------------##
## cd paths - disable here, let user choose in ~/.bashrc
##---------------------------------------------------------------------------##

CDPATH=

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
alias du='du -s -h'
alias free='free -m'
alias hosts='cat /etc/hosts'
alias hpss='echo Try using psi instead of hpss'
alias ldmesg='dmesg -s 65536 | less'
alias less='/usr/bin/less -r'
alias mdstat='cat /proc/mdstat'
alias meminfo='cat /proc/meminfo'
alias mroe='more'
alias print='lp'
alias resettermtitle='echo -ne "\033]0;${USER}@${target}: ${PWD}\007"'
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

# Use 'nt' to create new terminal window.
if test -z "$prefered_term"; then
   export prefered_term="gnome-terminal konsole xterm"
fi
if test "${TERM}" != emacs && 
   test "${TERM}" != dumb; then
   # choose a terminal {konsole,gnome-terminal,xterm}
   if test -z "$term"; then
     for item in ${prefered_term}; do
       tmp=`which ${item} 2> /dev/null`
       if test ${tmp:-notset} != notset; then
         export term=$item
         case ${item} in
         konsole)
           export term_opts='--vt_sz 90x60 --nomenubar'
           export title_flag=-T
           export exe_flag=-e
           break
           ;;
         gnome-terminal)
           export term_opts='--use-factory --geometry=80x60 --hide-menubar'
           export title_flag=-t
           export exe_flag=-x
           break
           ;;
         xterm)
           export term_opts='-geometry 80x60'
           export title_flag=-T
           export exe_flag=-e
           break
           ;;
         *)
           echo "You must must setup the variables term, term_opts,"
           echo "title_flag and exe_flag for this terminal ($item)"
           echo "in draco/environment/bashrc/.bashrc."  
           break
           ;;
         esac
       fi
     done
   fi
fi

# Function to create a new terminal window (requires X11)
if test ${term:-notset} != notset; then

function nt()
{
  ssh_cmd="ssh -AX"
  localmachine=`uname -n | sed -e 's/[.].*//'`
  machine=$1
  shift
  if test "${machine}notset" = notset ||
     test "${localmachine}" = "${machine}"; then
    machine=localhost
  fi
  case ${machine} in
  -h)
    opts="${title_flag} [machineName] ${exe_flag} ${ssh_cmd} [machineName]"
    echo "Usage: nt [machinename]"
    echo "  A new terminal window will be created using the command"
    echo "  ${term} ${term_opts} [extra args] $opts"
    return
    ;;
  localhost)
    machineName=localhost
    ;;
  *)
    machineName=${machine}
    ;;
  esac

  cmd="${term} ${term_opts} $* ${title_flag}"
  if test ${machineName} = localhost; then
    cmd="${cmd} ${localmachine}"
  else
    cmd="${cmd} ${machineName} ${exe_flag} ${ssh_cmd} ${machine}"
  fi
  echo $cmd
  eval $cmd
}
export nt
fi

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
## Machine specific settings
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

target="`uname -n | sed -e s/[.].*//`"
arch=`uname -m`

if test -d /ccs/codes/radtran; then # this is a CCS lan machine
   if test `uname -m` = "x86_64"; then
      target=ccslan64
   else
      target=ccslan32
   fi
fi

case $target in
ccslan32)
   source ${DRACO_SRC_DIR}/environment/bashrc/.bashrc_linux
   ;;
ccslan64)
   source ${DRACO_SRC_DIR}/environment/bashrc/.bashrc_linux64
   ;;
# flash / lightning
flash[a-d] | ffe-64 | lb-[1-7] | lc-64 | ll-[2-6] | lc-[1-6] | ffe[1-2] |\
    lb-dev)
    source ${DRACO_SRC_DIR}/environment/bashrc/.bashrc_bproc64
    ;;

# RoadRunner machines
rt-fe[1-4] | yr-fe1 | rt*[0-9]* | yra[0-9]* | rra[0-9]*a)
    source ${DRACO_SRC_DIR}/environment/bashrc/.bashrc_rr
    ;;
rr-dev-fe)
    source ${DRACO_SRC_DIR}/environment/bashrc/.bashrc_rr_dev
    ;;

#TLCC machines
tu-fe1 | tua[0-9]* | hu-fe[1-2] | hu*[0-9]*)
    source ${DRACO_SRC_DIR}/environment/bashrc/.bashrc_tlcc
    ;;

# Cielito
ct-fe1)
    source ${DRACO_SRC_DIR}/environment/bashrc/.bashrc_ct
    ;;

esac

source ${DRACO_SRC_DIR}/environment/bin/bash_functions.sh

##---------------------------------------------------------------------------##
## Aliases for machines
##---------------------------------------------------------------------------##

alias ccscs1='nt ccscs1'
alias ccscs8='nt ccscs8'
alias ccscs9='nt ccscs9'
alias rayo='nt rayo.lanl.gov'
alias rr='nt rr-dev-fe'
alias yr='nt yr-fe1'
alias turing='nt tu-fe1'

# No need to use ssh to pop a terminal from the current machine
alias ${target}='${term} ${term_opts}'

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

#Alternate
#PROMPT_COMMAND='DIR=`pwd|sed -e "s!$HOME!~!"`; if [ ${#DIR} -gt 30 ]; then CurDir=${DIR:0:12}...${DIR:${#DIR}-15}; else CurDir=$DIR; fi'
#PS1="[\$CurDir] \$ "

##---------------------------------------------------------------------------##
## Ensure that we have an ssh-agent running
##---------------------------------------------------------------------------##

 # if test -x /usr/bin/win-ssh-askpass.exe; then
 #   export SSH_ASKPASS=/usr/bin/win-ssh-askpass.exe
 # fi
 # if test -f ${HOME}/env.log; then
 #   rm -f ${HOME}/env.log
 # fi
 # set | grep SSH >& ${HOME}/env.log
 # if test -z ${SSH_AUTH_SOCK}; then
 #   if test -n "${verbose}"; then
 #      echo "no agent" 
 #   fi
 # else
 #   if test -n "`ssh-add -L | grep 'no identities'`"; then
 #      ssh-add < /dev/null
 #   fi
 # fi


##---------------------------------------------------------------------------##
## end of .bashrc
##---------------------------------------------------------------------------##
