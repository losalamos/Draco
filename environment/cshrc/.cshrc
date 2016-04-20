#!/bin/csh

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

# Extra module stuff
switch ("`uname -n`")
case pi-fey*.lanl.gov:
case pi*.localdomain:
case wf-fey*.lanl.gov
case wf*.localdomain:
case mu-fey*.lanl.gov:
case mu*.localdomain:
case mp*.lanl.gov:
case mp*.localdomain:
case ml-fey*.lanl.gov:
case ml*.localdomain:
    setenv VENDOR_DIR /usr/projects/draco/vendors
    module use $VENDOR_DIR/Modules/hpc
    module load friendly-testing user_contrib
###    module unload intel openmpi
    module load intel/15.0.5 openmpi
    module load git svn/1.8.10 emacs
    module load cmake/3.5.2 numdiff lapack/3.5.0 random123 eospac/6.2.4
    module load trilinos/12.6.1 superlu-dist/4.3
    module load parmetis/4.0.3 metis/5.1.0 ndi fstools totalview
    alias  topsn '/usr/projects/data/bin/latest/moonlight/topsn'
    breaksw

case lu*.lanl.gov:
case lu*.localdomain:
case ty*.lanl.gov:
case ty*.localdomain:
    setenv VENDOR_DIR /usr/projects/draco/vendors
    module use $VENDOR_DIR/Modules/hpc
    module use $VENDOR_DIR/Modules/tu-fe
    module load friendly-testing user_contrib
    module unload intel openmpi
    module load intel/15.0.3 openmpi
    module load git svn emacs
    module load cmake/3.5.2 numdiff lapack/3.5.0 random123 eospac/6.2.4
    module load trilinos/12.6.1 superlu-dist/4.3
    module load parmetis/4.0.3 ndi metis/5.1.0
    alias  topsn '/usr/projects/data/bin/latest/moonlight/topsn'
    breaksw

case redfta[0-9]*:
case rfta*:
    setenv VENDOR_DIR /usr/projects/draco/vendors
    module use $VENDOR_DIR/Modules/hpc
    module use $VENDOR_DIR/Modules/tu-fe
    module load git svn
    breaksw

case ct*:
case ci*:
    setenv VENDOR_DIR /usr/projects/draco/vendors
    # source /usr/projects/crestone/dotfiles/Cshrc
    module use $VENDOR_DIR/Modules/hpc
    module load user_contrib friendly-testing

    # Move some environment out of the way.
    module unload PrgEnv-intel PrgEnv-pgi
    module unload cmake numdiff svn gsl
    module unload papi perftools

    # load the Intel programming env, but then unloda libsci and totalview
    module load PrgEnv-intel # this loads xt-libsci and intel/XXX
    module unload xt-libsci intel # xt-totalview
    module load intel/14.0.4.211

    # draco modules start here.
    module load gsl/1.15 lapack/3.5.0
    module load cmake/3.3.1 numdiff svn git emacs
    module load trilinos SuperLU_DIST
    module load ParMetis ndi random123 eospac/6.2.4

    setenv OMP_NUM_THREADS 8

    # Avoid run time messages of the form:
    # "OMP: Warning #72: KMP_AFFINITY: affinity only supported for Intel(R) processors."
    # Ref: http://software.intel.com/en-us/articles/bogus-openmp-kmp_affinity-warnings-on-non-intel-processor-hosts/
    setenv KMP_AFFINITY none
    breaksw

case seqlac*:
    setenv VENDOR_DIR /usr/gapps/jayenne/vendors

    # LLNL uses dotkit instead of modules
    setenv DK_NODE ${DK_NODE}:${VENDOR_DIR}/Modules/sq

    # Draco dotkits
    use xlc12 # use gcc472
    use numdiff
    use random123

    # LLNL dotkits
    use cmake331
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
