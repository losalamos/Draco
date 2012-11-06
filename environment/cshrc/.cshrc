#!/bin/csh

# Use: In ~/.cshrc add the following code:
#
# setenv DRACO_ENV_DIR ~/draco/environment
# source $DRACO_ENV_DIR/cshrc/.cshrc
#

setenv PATH $DRACO_ENV_DIR/bin:$PATH

# Extra module stuff
switch ("`uname -n`")
case pi-fey*.lanl.gov:
case pi*.localdomain:

#    source /usr/projects/crestone/dotfiles/Cshrc
    source /usr/projects/draco/vendors/modules-3.2.9/init/tcsh

    module use $DRACO_ENV_DIR/Modules/hpc
    module use $DRACO_ENV_DIR/Modules/tu-fe
    module load friendly-testing 
    module load intel/12.1.5 openmpi
    module load gsl/1.14-intel svn emacs
    module load cmake numdiff git lapack/3.4.1-intel
    module load trilinos SuperLU_DIST
    module load ParMetis ndi random123 eospac
    # PGI keeps running out of tmp sapce
#     setenv TMPDIR /scratch/$USER/tmp
#     if (! -d $TMPDIR ) then
#        mkdir $TMPDIR
#     endif
    breaksw
case lu*.lanl.gov:
case lu*.localdomain:
case ty*.lanl.gov:
case ty*.localdomain:
    module use $DRACO_ENV_DIR/Modules/hpc
    module use $DRACO_ENV_DIR/Modules/tu-fe
    module load friendly-testing 
    module load intel/12.1.5 openmpi
    module load gsl/1.14-intel svn emacs
    module load cmake numdiff git lapack/3.4.1-intel
    module load trilinos SuperLU_DIST
    module load ParMetis ndi
    breaksw

case ml-fey*.lanl.gov:
case ml*.localdomain:

    if ($?tcsh) then
        set modules_shell="tcsh"
    else
        set modules_shell="csh"
    endif
    
    set exec_prefix='/usr/projects/draco/vendors/modules-3.2.9'
    set prefix=""
    set postfix=""

    if ( $?histchars ) then
        set histchar = `echo $histchars | cut -c1`
        set _histchars = $histchars
        set prefix  = 'unset histchars;'
        set postfix = 'set histchars = $_histchars;'
    else
        set histchar = \!
    endif
    
    if ($?prompt) then
        set prefix  = "$prefix"'set _prompt="$prompt";set prompt="";'
        set postfix = "$postfix"'set prompt="$_prompt";unset _prompt;'
    endif
    
    if ($?noglob) then
        set prefix  = "$prefix""set noglob;"
        set postfix = "$postfix""unset noglob;"
    endif
    set postfix = "set _exit="'$status'"; $postfix; test 0 = "'$_exit;'
    
    alias module $prefix'eval `'$exec_prefix'/bin/modulecmd '$modules_shell' '$histchar'*`; '$postfix
    unset exec_prefix
    unset prefix
    unset postfix
    
    setenv MODULESHOME /usr/projects/draco/vendors/modules-3.2.9
    setenv MODULE_VERSION 3.2.9
    setenv MODULE_VERSION_STACK 3.2.9

    source /usr/projects/draco/vendors/modules-3.2.9/init/tcsh

    module use $DRACO_ENV_DIR/Modules/hpc
    module use $DRACO_ENV_DIR/Modules/tu-fe
    module load friendly-testing 
    module load intel/12.1.5 openmpi cudatoolkit
    module load cmake gsl/1.14-intel svn 
    module load numdiff lapack/3.4.1-intel
    module load trilinos SuperLU_DIST/3.0-intel
    module load ParMetis/3.1.1-intel ndi random123 eospac
    alias  mvcap 'cd /usr/projects/capsaicin/devs/jhchang'  
    breaksw
case ct*:
case ci*:
   # source /usr/projects/crestone/dotfiles/Cshrc
   module use $DRACO_ENV_DIR/Modules/hpc
   module use $DRACO_ENV_DIR/Modules/ct-fe
   # Move some environment out of the way.
   module unload PrgEnv-intel PrgEnv-pgi
   module unload cmake numdiff svn gsl
   # load the Intel programming env, but then unloda libsci and totalview
   module load PrgEnv-intel
   module unload xt-libsci xt-totalview 
   # draco modules start here.
   module load gsl/1.14 lapack/3.4.1-intel
   module load cmake numdiff svn
   module load trilinos SuperLU_DIST/3.0-intel 
   module load ParMetis/3.1.1-intel ndi random123 eospac

   # Avoid run time messages of the form:
   # "OMP: Warning #72: KMP_AFFINITY: affinity only supported for Intel(R) processors."
   # Ref: http://software.intel.com/en-us/articles/bogus-openmp-kmp_affinity-warnings-on-non-intel-processor-hosts/
   setenv KMP_AFFINITY none
   breaksw
case rr-dev*:
case rra[0-9][0-9][0-9]a*:
   source /usr/projects/eap/dotfiles/Cshrc
   module use $DRACO_ENV_DIR/Modules/hpc
   module use $DRACO_ENV_DIR/Modules/rr-dev-fe
   module load friendly-testing cellsdk svn
   module unload pgi openmpi-pgi
   module load cmake numdiff python openmpi-gcc/1.4.3
   module load gcc/4.7.1
   breaksw
case rra[0-9][0-9][0-9][bcd]*:
   # source /usr/projects/crestone/dotfiles/Cshrc
   module use $DRACO_ENV_DIR/Modules/ppc64
   module load friendly-testing cellsdk
   module load cmake gsl-1.14 numdiff 
   module load 
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
