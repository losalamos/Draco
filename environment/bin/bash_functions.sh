## -*- Mode: sh -*-
##---------------------------------------------------------------------------##
## File  : environment/bin/bash_functions.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##
##
## Summary: Misc bash functions useful during development of code.
##
## Functions
## ---------
##
## whichall <string> - show path of all commands in path that match
##                     <string>
##
## cleanemacs        - recursively remove ~ files, .flc files and .rel
##                     files.
##
## findsymbol <sym>  - search all libraries (.so and .a files) in the
##                     current directory for symbol <sym>.
##
## npwd              - function used to set the prompt under bash.
##
## xfstatus          - report status of transfer.lanl.gov
##
## rm_from_path      - remove a directory from $PATH
##
## add_to_path       - add a directory to $PATH
##
## proxy             - (un)set http_proxy variables
##
## fn_exists         - test if a bash function is defined
##
## run               - echo then evaluate a bash command
##
## rdde              - reload the default draco environment
##
## qrm               - quick remove (for lustre filesystems).
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## Find all matches in PATH (not just the first one)
##---------------------------------------------------------------------------##

function whichall ()
{
  for dir in ${PATH//:/ }; do
    if [ -x $dir/$1 ]; then
      echo $dir/$1;
    fi;
  done
}

##---------------------------------------------------------------------------##
## Recursively delete all ~ files.
##---------------------------------------------------------------------------##

function cleanemacs
{
  echo "Cleaning up XEmacs temporary and backup files.";
  find . -name '*~' -exec echo rm -rf {} \;
  find . -name '*~' -exec rm -rf {} \;
  find . -name '.*~' -exec echo rm -rf {} \;
  find . -name '.*~' -exec rm -rf {} \;
  find . -name '*.flc' -exec echo rm -rf {} \;
  find . -name '*.flc' -exec rm -rf {} \;
  find . -name '*.rel' -exec echo rm -rf {} \;
  find . -name '*.rel' -exec rm -rf {} \;
  if test -d "${HOME}/.emacs-flc"; then
    echo "rm -r ${HOME}/.emacs-flc";
    rm -r ${HOME}/.emacs-flc;
  fi;
  echo "done."
}

##---------------------------------------------------------------------------##
## Used for formatting PROMPT.
## $HOME            -> ~
## ...scratch...    -> #
## .../projects/... -> @
##---------------------------------------------------------------------------##
parse_git_branch() {
  git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/ (\1)/'
}

function npwd()
{
  # Optional arguments:
  #   $1 - number of chars to print.
  #   $2 - scratch location

  # How many characters of the $PWD should be kept
  local pwdmaxlen=40
  if [[ $1 ]]; then pwdmaxlen=$1; fi

  # local regHome=$(echo ${HOME} | sed -e 's/.*\///')

  local scratchdirs=/scratch:/lustre/ttscratch1:/lustre/scratch[123]/yellow
  if [[ $2 ]]; then scratchdirs=$2; fi

  # Indicator that there has been directory truncation:
  local trunc_symbol="..."
  # substitute ~ for $HOME to shorten the full path
  newPWD=$(echo ${PWD} | sed -e "s%$HOME%~%")
  local oldIFS=$IFS
  IFS=:
  for dir in $scratchdirs; do
    newPWD=$(echo ${newPWD} | sed -e "s%${dir}/${USER}%#%")
  done
  IFS=$oldIFS

  local devdirs=/usr/projects/jayenne/devs
  newPWD=$(echo ${newPWD} | sed -e "s%${devdirs}/${USER}%@%")

  if [ ${#newPWD} -gt $pwdmaxlen ]; then
    local pwdoffset=$(( ${#newPWD} - $pwdmaxlen ))
    newPWD="${trunc_symbol}${newPWD:$pwdoffset:$pwdmaxlen}"
  fi
  echo ${newPWD}
}

# Alternate
function npwd_alt()
{
  # Home directory
  local regHome=${HOME}
  # Home directory might be prepended with NFS mount location.
  local netHome=/net/machname.lanl.gov
  #   How many characters of the $PWD should be kept
  local pwdmaxlen=40
  #   Indicator that there has been directory truncation:
  local trunc_symbol="..."
  # substitute ~ for $HOME to shorten the full path
  newPWD=`echo ${PWD} | sed -e "s%${netHome}%%" | sed -e "s%${regHome}%~%"`
  if test ${#newPWD} -gt $pwdmaxlen ; then
    local pwdoffset=$(( ${#newPWD} - $pwdmaxlen ))
    newPWD="${trunc_symbol}${newPWD:$pwdoffset:$pwdmaxlen}"
  fi
  echo ${newPWD}
}
# Another alternate
# PROMPT_COMMAND='DIR=`pwd|sed -e "s!$HOME!~!"`; if [ ${#DIR} -gt 30 ]; then CurDir=${DIR:0:12}...${DIR:${#DIR}-15}; else CurDir=$DIR; fi'
# PS1="[\$CurDir] \$ "

##---------------------------------------------------------------------------##
## Usage:
##    findsymbol <symbol>
##
## Searches all .a and .so files in local directory for symbol
## <symbol>.  If found, the script provides the name of the library
## that contains the symbol.
##---------------------------------------------------------------------------##

function findsymbol()
{
  local nm_opt='-a'
  local a_libs=`\ls -1 *.a`
  if test -z "$a_libs"; then a_libs=""; fi
  local so_libs=`\ls -1 *.so`
  if test -z "$so_libs"; then so_libs=""; fi
  local libs="$a_libs $so_libs"
  echo " "
  echo "Searching..."
  local symbol=" T "
  for lib in $libs; do
    local gres=`nm $nm_opt $lib | grep $1 | grep "$symbol"`
    if ! test "$gres" = ""; then
      echo " "
      echo "Found \"$symbol\" in $lib:"
      echo "     $gres"
    fi
  done
  echo " "
}

##---------------------------------------------------------------------------##
## Transfer 2.0 (Mercury replacement)
## Ref: http://transfer.lanl.gov
##
## Examples:
##   xfpush foo.txt
##   xfstatus
##   xfpull foo.txt
##---------------------------------------------------------------------------##

function xfstatus()
{
  ssh red@transfer.lanl.gov myfiles
}

##---------------------------------------------------------------------------##
## If string is found in PATH, remove it.
##---------------------------------------------------------------------------##
function rm_from_path ()
{
  badpath=$1
  newpath=""
  for dir in ${PATH//:/ }; do
    if ! test "${badpath}" = "${dir}"; then
      newpath="${newpath}:${dir}"
    fi;
  done
  newpath=`echo $newpath | sed -e s/^[:]//`
  export PATH=$newpath
}

##---------------------------------------------------------------------------##
## If path is a directory add it to PATH (if not already in PATH)
##
## Use:
##   add_to_path <path> TEXINPUTS|BSTINPUTS|BIBINPUTS|PATH
##---------------------------------------------------------------------------##
function add_to_path ()
{
  case $2 in
    TEXINPUTS)
      if [ -d "$1" ] && [[ ":${TEXINPUTS}:" != *":$1:"* ]]; then
        TEXINPUTS="${TEXINPUTS:+${TEXINPUTS}:}$1"; fi ;;
    BSTINPUTS)
      if [ -d "$1" ] && [[ ":${BSTINPUTS}:" != *":$1:"* ]]; then
        BSTINPUTS="${BSTINPUTS:+${BSTINPUTS}:}$1"; fi ;;
    BIBINPUTS)
      if [ -d "$1" ] && [[ ":${BIBINPUTS}:" != *":$1:"* ]]; then
        BIBINPUTS="${BIBINPUTS:+${BIBINPUTS}:}$1"; fi ;;
    *)
      if [ -d "$1" ] && [[ ":${PATH}:" != *":$1:"* ]]; then
        PATH="${PATH:+${PATH}:}$1"; fi ;;
  esac
}

##---------------------------------------------------------------------------##
## Toggle LANL proxies on/off
## https://wiki.archlinux.org/index.php/proxy_settings
##---------------------------------------------------------------------------##
function proxy()
{
  if [[ ! ${http_proxy} ]]; then
    # proxies not set, set them
    export http_proxy=http://proxyout.lanl.gov:8080
    export https_proxy=$http_proxy
    export HTTP_PROXY=$http_proxy
    export HTTPS_PROXY=$http_proxy
    # export http_no_proxy="*.lanl.gov"
    export no_proxy="localhost,127.0.0.1,.lanl.gov"
    export NO_PROXY=$no_proxy
  else
    # proxies are set, kill them
    unset http_proxy
    unset https_proxy
    unset HTTP_PROXY
    unset HTTPS_PROXY
    #unset http_no_proxy
    unset no_proxy
    unset NO_PROXY
  fi
}

##---------------------------------------------------------------------------##
## Test to determine if named bash function exists in the current environment.
##---------------------------------------------------------------------------##
function fn_exists()
{
  type $1 2>/dev/null | grep -q 'is a function'
  res=$?
  echo $res
  return $res
}

##---------------------------------------------------------------------------##
## Echo commands before execution (used in scripts)
##---------------------------------------------------------------------------##
function run () {
  echo $1
  if ! [ $dry_run ]; then eval $1; fi
}

##---------------------------------------------------------------------------##
## Reset the draco developer environment
##---------------------------------------------------------------------------##
function rdde ()
{
  unset DRACO_BASHRC_DONE
  source ${DRACO_ENV_DIR}/bashrc/.bashrc
}

#------------------------------------------------------------------------------#
# Quick remove: instead of 'rm -rf', mv the directory to .../trash/tmpname
#
# Use: qrm dir1 dir2 ...
#------------------------------------------------------------------------------#
function qrm ()
{
  # must provide at least one directory
  if [[ ! $1 ]]; then
    echo "You must provide a single argument to this function that is the name of the directory to delete."
    return
  fi

  for dir in "$@"; do

    # fully qualified directory name
    fqd=`cd $dir && pwd`

    if [[ ! -d $dir ]]; then
      echo "$dir doesn't look like a directory. aborting."
      return
    fi

    if [[ `echo $fqd | grep -c scratch` == 0 ]]; then
      echo "This command should only be used for scratch directories."
      return
    fi

    # Identify the scratch system
    trashdir=`echo $fqd | sed -e "s%$USER.*%$USER/trash%"`

    # ensure trash folder exists.
    mkdir -p $trashdir
    if ! [[ -d $trashdir ]]; then
      echo "FATAL ERROR: Unable access trashdir = $trashdir"
      return
    fi

    # We rename/move the old directory to a random name
    TMPDIR=$(mktemp -d $trashdir/XXXXXXXXXX) || { echo \
      "Failed to create temporary directory"; return; }

    # do it
    mv $dir $TMPDIR/.

  done
}

#------------------------------------------------------------------------------#
# End environment/bin/.bash_functions
#------------------------------------------------------------------------------#
