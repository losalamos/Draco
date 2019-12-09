#!/bin/csh
##---------------------------------------------------------------------------##
## File  : environment/cshrc/.cshrc
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2017, Triad National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# Use: In ~/.cshrc add the following code:
#
# setenv DRACO_ENV_DIR ~/draco/environment
# source $DRACO_ENV_DIR/cshrc/.cshrc
#

if( ! $?DRACO_ENV_DIR )then
  if( -d /usr/projects/draco/environment )then
    setenv DRACO_ENV_DIR /usr/projects/draco/environment
  endif
  if( -d ~/draco/environment )then
    setenv DRACO_ENV_DIR ~/draco/environment
  endif
endif
if( ! $?DRACO_ENV_DIR )then
  echo "ERROR: DRACO_ENV_DIR has not been set."
  exit 1
endif
if( ! -d "${DRACO_ENV_DIR}" )then
  echo "ERROR: DRACO_ENV_DIR is not a valid directory (${DRACO_ENV_DIR})"
  exit 1
endif

# Setup PATH
setenv PATH $DRACO_ENV_DIR/bin:$PATH
if( -d ~/bin )then
    setenv PATH ~/bin:$PATH
endif

# Hooks for clang-format as git commit hook:
# Possible values: ON, TRUE, OFF, FALSE, DIFF (the default value is ON).
setenv DRACO_AUTO_CLANG_FORMAT ON

# Extra module stuff
switch ("`uname -n`")
case pi-fey*.lanl.gov:
case pi*.localdomain:
case wf-fey*.lanl.gov
case wf*.localdomain:
case sn*:
case fi*:
case ic*:
    setenv VENDOR_DIR /usr/projects/draco/vendors
    module load friendly-testing user_contrib
    module load clang-format
    module load intel/17.0.1 openmpi/1.10.5 mkl
    module load git subversion emacs grace
    module load cmake/3.6.2 numdiff random123 eospac/6.2.4
    module load trilinos/12.8.1 superlu-dist/4.3
    module load parmetis/4.0.3 metis/5.1.0 ndi totalview
    alias  topsn '/usr/projects/data/bin/latest/moonlight/topsn'
    breaksw

case redfta[0-9]*:
case rfta*:
case redcap*:
    setenv VENDOR_DIR /usr/projects/draco/vendors
    module use $VENDOR_DIR/Modules/tu-fe
    module load git svn
    breaksw

case tt*:
case tr*:

    setenv VENDOR_DIR /usr/projects/draco/vendors
    # source /usr/projects/crestone/dotfiles/Cshrc
    module load user_contrib friendly-testing

    # Move some environment out of the way.
    module unload PrgEnv-intel PrgEnv-pgi
    module unload cmake numdiff svn gsl
    module unload papi perftools
    module load clang-format

    # load the Intel programming env, but then unloda libsci and totalview
    module load PrgEnv-intel # this loads xt-libsci and intel/XXX
    module unload intel
    module load intel/17.0.1
    module unload cray-libsci gcc/6.1.0

    # draco modules start here.
    module load metis parmetis/4.0.3 trilinos/12.8.1 superlu-dist/4.3
    module load gsl/2.1 cmake/3.6.2 numdiff ndi random123 eospac/6.2.4
    module load subversion git

    setenv OMP_NUM_THREADS 16
    setenv CXX CC
    setenv CC cc
    setenv FC ftn
    setenv CRAYPE_LINK_TYPE dynamic

    breaksw

case seqlac*:
    setenv VENDOR_DIR /usr/gapps/jayenne/vendors

    # LLNL uses dotkit instead of modules
    setenv DK_NODE ${DK_NODE}:${VENDOR_DIR}/Modules/sq

    # Draco dotkits
    # use xlc12
    use gcc484
    use numdiff
    use random123
    use gsl

    # LLNL dotkits
    use cmake361
    use erase=del
    use alia1++

    unalias rm

    breaksw

endsw

# Set term title
set myhost=`echo $HOST | sed -e 's/[.].*//g'`

# Aliases
alias btar 'tar --use-compress-program /usr/bin/bzip2'
alias cd.. 'cd ..'
alias cpuinfo 'cat /proc/cpuinfo'
alias df 'df -h'
alias dirs 'dirs -v'
alias dmesg 'dmesg -s 65536'
alias du 'du -s -h'
alias em 'emacs $* -g 90x55'
alias free 'free -m'
alias hosts 'cat /etc/hosts'
alias hpss 'echo Try using psi instead of hpss'
alias l. 'ls --color -aFl'
alias ldmesg 'dmesg -s 65536 | less'
alias less '/usr/bin/less -r'
alias ll 'ls --color -Fl'
alias ls 'ls --color -F'
alias lt 'ls --color -Flt'
alias lt. 'ls --color -aFlt'
alias mdstat 'cat /proc/mdstat'
alias meminfo 'cat /proc/meminfo'
alias mroe 'more'
alias print 'lp'
alias resettermtitle 'echo -ne "\033]0;${myhost}\007"'
alias rtt 'echo -ne "\033]0;${myhost}\007"'
alias sdk 'export DISPLAY 128.165.87.170:0.0;echo The DISPLAY value is now: $DISPLAY'
alias sdl 'DISPLAY 127.0.0.1:0.0;echo The DISPLAY value is now: $DISPLAY'
alias vi 'vim'
alias watchioblocks 'ps -eo stat,pid,user,command | egrep "^STAT|^D|^R"'
alias which 'alias | /usr/bin/which --tty-only --read-alias --show-dot --show-tilde'
alias wmdstat 'watch -n 2 "cat /proc/mdstat"'
alias xload 'xload -label ${myhost}'
