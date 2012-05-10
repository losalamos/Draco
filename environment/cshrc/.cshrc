# -*- mode:sh -*-

# Use: In ~/.cshrc add the following code:
#
# setenv DRACO_ENV_DIR ~/draco/environment
# source $DRACO_ENV_DIR/cshrc/.cshrc
#

setenv PATH $DRACO_ENV_DIR/bin:$PATH

# Extra module stuff
switch ("`uname -n`")
case tu*.lanl.gov:
case tu*.localdomain:
case hu*.lanl.gov:
case hu*.localdomain:
case ty*.lanl.gov:
case ty*.localdomain:
    source /usr/projects/crestone/dotfiles/Cshrc
    module use $DRACO_ENV_DIR/Modules/hpc
    module use $DRACO_ENV_DIR/Modules/tu-fe
    module load friendly-testing 
    module load intel-c intel-f openmpi-intel
    module load gsl/1.14-intel svn emacs
    module load cmake numdiff git xshow papi/4.1.3
    # PGI keeps running out of tmp sapce
    setenv TMPDIR /scratch/$USER/tmp
    if (! -d $TMPDIR ) then
       mkdir $TMPDIR
    endif
    breaksw
case yr*.lanl.gov:
case yr*:
    source /usr/projects/crestone/dotfiles/Cshrc
    module use $DRACO_ENV_DIR/Modules/hpc
    module use $DRACO_ENV_DIR/Modules/yr-fe
    module load lapack/atlas-3.8.3 svn
    module load cmake numdiff git xshow
    module load gsl/1.14-pgi emacs
    # PGI keeps running out of tmp sapce
    setenv TMPDIR /scratch/$USER/tmp
    if (! -d $TMPDIR ) then
       mkdir $TMPDIR
    endif
    breaksw
case ct*:
case ci*:
   # source /usr/projects/crestone/dotfiles/Cshrc
   module use $DRACO_ENV_DIR/Modules/hpc
   module use $DRACO_ENV_DIR/Modules/ct-fe
   module load gsl/1.14 svn
   module load cmake/2.8.7 numdiff git xshow papi
   module load tkdiff/4.1.4 openspeedshop/2.0.1b10 
   module unload xt-libsci
   # module load lapack/3.4.0-pgi # use /opt/pgi/11.10.0/...
   breaksw
case rr-dev*:
case rra[0-9][0-9][0-9]a*:
   source /usr/projects/crestone/dotfiles/Cshrc
   module use $DRACO_ENV_DIR/Modules/hpc
   module use $DRACO_ENV_DIR/Modules/rr-dev-fe
   module load friendly-testing cellsdk svn
   module unload pgi openmpi-pgi
   module load cmake numdiff git xshow python openmpi-gcc/1.4.3
   breaksw
case rra[0-9][0-9][0-9][bcd]*:
   # source /usr/projects/crestone/dotfiles/Cshrc
   module use $DRACO_ENV_DIR/Modules/ppc64
   module load friendly-testing cellsdk
   module load cmake gsl-1.14 numdiff 
   module load 
   breaksw
case lu*.lanl.gov
case lu*.localdomain:
    source /usr/projects/crestone/dotfiles/Cshrc
    module use $DRACO_ENV_DIR/Modules/hpc
    module use $DRACO_ENV_DIR/Modules/tu-fe
    module load friendly-testing 
    module load intel-c intel-f openmpi-intel
    module load gsl/1.14-intel svn emacs
    module load cmake numdiff
    breaksw
case gondolin*:
    source /ccs/codes/radtran/vendors/modules-3.2.7/init/csh
    module load grace totalview numdiff git gsl svn gcc lapack/3.4.0
    module load cmake openmpi emacs
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
