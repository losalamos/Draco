#!/bin/csh

# Use: In ~/.cshrc add the following code:
#
 setenv DRACO_ENV_DIR ~/draco/environment
# source $DRACO_ENV_DIR/cshrc/.cshrc
#

setenv PATH $DRACO_ENV_DIR/bin:$PATH

# Extra module stuff
switch ("`uname -n`")
case pi-fey*.lanl.gov:
case pi*.localdomain:
case mu-fey*.lanl.gov:
case mu*.localdomain:

#    source /usr/projects/crestone/dotfiles/Cshrc
    source /usr/projects/draco/vendors/modules-3.2.9/init/tcsh

    module use $DRACO_ENV_DIR/Modules/hpc
    module use $DRACO_ENV_DIR/Modules/tu-fe
    module load friendly-testing 
    module load intel/13.1.0 openmpi/1.6.3
    module load gsl/1.14-intel svn emacs
    module load cmake numdiff git lapack/3.4.1-intel
    module load trilinos SuperLU_DIST
    module load ParMetis ndi random123 eospac
    # PGI keeps running out of tmp sapce
#     setenv TMPDIR /scratch/$USER/tmp
#     if (! -d $TMPDIR ) then
#        mkdir $TMPDIR
#     endif
    setenv VENDOR_DIR /usr/projects/draco/vendors
    breaksw
case lu*.lanl.gov:
case lu*.localdomain:
case ty*.lanl.gov:
case ty*.localdomain:
    module unload openmpi-intel
    module use $DRACO_ENV_DIR/Modules/hpc
    module use $DRACO_ENV_DIR/Modules/tu-fe
    module load friendly-testing 
    module load intel/13.1.0 openmpi/1.6.3
    module load gsl/1.14-intel emacs
    module load cmake numdiff git lapack/3.4.1-intel
    module load trilinos SuperLU_DIST
    module load ParMetis ndi
    setenv VENDOR_DIR /usr/projects/draco/vendors
    breaksw

case redfta[0-9]*:
    module use $DRACO_ENV_DIR/Modules/hpc
    module use $DRACO_ENV_DIR/Modules/tu-fe
    module load git svn
    breaksw

case ml-fey*.lanl.gov:
case ml*.localdomain:
    module use $DRACO_ENV_DIR/Modules/hpc
    module use $DRACO_ENV_DIR/Modules/tu-fe
    module load friendly-testing 
    module load intel/13.1.0 openmpi/1.6.3 # cudatoolkit
    module load cmake gsl/1.14-intel svn fstools 
    module load numdiff lapack/3.4.1-intel totalview
    module load SuperLU_DIST/3.0-openmpi163-intel1310
    module load trilinos/10.12.2-openmpi163-intel1310
    module load ParMetis/3.1.1-openmpi163-intel1310 
    module load ndi random123 eospac/v6.2.4beta.1-moonlight
    alias  topsn '/usr/projects/data/bin/latest/moonlight/topsn' 
    setenv VENDOR_DIR /usr/projects/draco/vendors
    breaksw

case ct*:
case ci*:
   # source /usr/projects/crestone/dotfiles/Cshrc
   module use $DRACO_ENV_DIR/Modules/hpc
   module use $DRACO_ENV_DIR/Modules/ct-fe
   # Move some environment out of the way.
   module unload PrgEnv-intel PrgEnv-pgi
   module unload cmake numdiff svn gsl
   module unload papi perftools
   # load the Intel programming env, but then unloda libsci and totalview
   module load PrgEnv-intel
   module unload xt-libsci xt-totalview intel
   module load intel/13.0.1.117
   # draco modules start here.
   module load gsl/1.14 lapack/3.4.1-intel
   module load cmake numdiff subversion emacs
   module load trilinos SuperLU_DIST/3.0-intel 
   module load ParMetis/3.1.1-intel ndi random123 eospac/v6.2.4beta.1-cielito

   setenv OMP_NUM_THREADS 8

   # Avoid run time messages of the form:
   # "OMP: Warning #72: KMP_AFFINITY: affinity only supported for Intel(R) processors."
   # Ref: http://software.intel.com/en-us/articles/bogus-openmp-kmp_affinity-warnings-on-non-intel-processor-hosts/
   setenv KMP_AFFINITY none
   setenv VENDOR_DIR /usr/projects/draco/vendors
   breaksw
case gondolin*:
    source /ccs/codes/radtran/vendors/modules-3.2.7/init/csh
    module load grace totalview numdiff git gsl svn gcc lapack/3.4.0
    module load cmake openmpi emacs random123
    module load trilinos BLACS SCALAPACK SuperLU_DIST hypre/2.0.0 ndi 
    module load ParMetis/3.1.1
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
