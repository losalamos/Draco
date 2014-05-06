## -*- Mode: sh -*-----------------------------------------------------------##
##
## file  : bash_functions.sh
##
## Summary: Misc bash functions useful during development of code.
##
## 1. Use GNU tools instead of vendor tools when possible
## 2. Create some alias commands to provide hints when invalid
##    commands are issued.
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
## ssh1/scp1         - force use of protocol version 1.
##
## findsymbol <sym>  - search all libraries (.so and .a files) in the
##                     current directory for symbol <sym>.
##
## pkgdepends        - Print a list of dependencies for the current
##                     directory. 
##
## npwd              - function used to set the prompt under bash.
##
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## EOSPAC Setup
##---------------------------------------------------------------------------##

if test -d /ccs/codes/radtran/physical_data/eos; then
   export SESPATHU=/ccs/codes/radtran/physical_data/eos
   export SESPATHC=/ccs/codes/radtran/physical_data/eos
elif test -d /usr/projects/data/eos; then
   export SESPATHU=/usr/projects/data/eos
   export SESPATHU=/usr/projects/data/eos
fi

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

function ssh1 { ssh -x -1 $*; }
function scp1 { scp -oProtocol=1 $*; }

##---------------------------------------------------------------------------##
## Used for formatting PROMPT.
##---------------------------------------------------------------------------##

function npwd() 
{
   local regHome=$(echo ${HOME} | sed -e 's/.*\///')
   #   How many characters of the $PWD should be kept
   local pwdmaxlen=40
   #   Indicator that there has been directory truncation:
   local trunc_symbol="..."
   # substitute ~ for $HOME to shorten the full path
   newPWD=$(echo ${PWD} | sed -e "s/.*${regHome}/~/")
   if [ ${#newPWD} -gt $pwdmaxlen ]
   then
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
    if test `uname` = OSF1; then
       nm_opt=''
    fi
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
## Usage:
##    pkgdepends
##
## Purpose:
##    The script will list all of the vendors and Draco packages
##    dependencies for the files in the current directory.
##---------------------------------------------------------------------------##

function pkgdepends()
{
   echo "This package depends on:"
   echo " "
   echo "Packages:"
   grep 'include [<"].*[/]' *.cc *.hh | sed -e 's/.*[#]include [<"]/   /' | sed -e 's/\/.*//' | sort -u
   echo " "
   if test -f configure.ac; then
      echo "Vendors:"
      grep "SETUP[(]pkg" configure.ac | sed -e 's/AC_/   /' | sed -e 's/_.*//'
   fi
}

##---------------------------------------------------------------------------##
## Usage:
##    findgrep <regex>
##
## Finds all occurances of <regex> by looking at all files in the
## current directory and all subdirectories recursively.  Prints the
## filename and line number for each occurance.
##---------------------------------------------------------------------------##
function findgrep()
{
    # Exclude .svn directories: (-path '*/.svn' -prune)
    # Or (-o)
    files=`find . -path '*/.svn' -prune -o -type f -exec grep -q $1 {} /dev/null \; -print`
    for file in $files; do
      echo " "
      echo "--> Found \"$1\" in file \"$file\" :"
      echo " "
      grep $1 $file
    done
    echo " "
}

##---------------------------------------------------------------------------##
## Usage:
##    archive [age_in_days]
##
## Move all files older than [age_in_days] (default: 7d) from the
## current directory into a subdirectory named as the current year.
##---------------------------------------------------------------------------##
function archive()
{
  # Find files (-type f) older than 7 days (-mtime +7) and in the local
  # directory (-maxdepth 1) and move them to a subdirectory named by
  # the current year.
  local dt=7
  if test -n "$1"; then
    dt=$1
  fi
  local year=`date +%Y`
#  local year=`stat {} | grep Change | sed -e 's/Change: //' -e 's/[-].*//'
#  if ! test -d $year; then
#    mkdir $year
#  fi
#  echo "Moving files to ${year}/..."
#  cmd="find . -maxdepth 1 -mtime +${dt} -type f -exec mv {} ${year}/. \;"
  # echo $cmd
#  eval $cmd
  files=`find . -maxdepth 1 -mtime +${dt} -type f`
  for file in $files; do
    year=`stat ${file} | grep Modify | sed -e 's/Modify: //' -e 's/[-].*//'`
    echo "   Moving $file to ${year}"
    if ! test -d $year; then
      mkdir $year
    fi
    mv ${file} $year/$file
  done
}

##---------------------------------------------------------------------------##
## Use 'wiki <term>'
##---------------------------------------------------------------------------##
function wiki()
{
  dig +short txt $1.wp.dg.cx;
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
function xfpush()
{
    myfiles="$*"
    if ! test -n "$1"; then
       echo "ERROR: You must profile a file for transfering: xfpush foo.txt"
       return
    fi
    for myfile in $myfiles; do
       scp $myfile red@transfer.lanl.gov:
       echo scp $myfile red@transfer.lanl.gov:
    done
}
function xfstatus()
{
    ssh red@transfer.lanl.gov myfiles
}
function xfpull()
{
    wantfiles="$*"
    filesavailable=`ssh red@transfer.lanl.gov myfiles`
    for wantfile in $wantfiles; do

    # sanity check: is the requested file in the list?
    fileready=`echo $filesavailable | grep $wantfile`
    if test "${fileready}x" = "x"; then
        echo "ERROR: File '${wantfile}' is not available (yet?) to pull."
        echo "       Run 'xfstatus' to see list of available files."
        return
    fi
    # Find the file identifier for the requested file.  The variable
    # filesavailable contains a list of pairs:  
    # { (id1, file1), (id2, file2), ... }.  Load each pair and if the
    # filename matches the requested filename then pull that file id.
    # Once pulled, remove the file from transfer.lanl.gov.
    is_file_id=1
    for entry in $filesavailable; do
        if test $is_file_id = 1; then
            fileid=$entry
            is_file_id=0
        else
            if test $entry = $wantfile; then
                echo "scp red@transfer.lanl.gov:${fileid} ."
                scp red@transfer.lanl.gov:${fileid} .
                echo "ssh red@transfer.lanl.gov delete ${fileid}"
                ssh red@transfer.lanl.gov delete ${fileid}
                return
            fi
            is_file_id=1
        fi
    done

    done # end loop over $wantfiles
}


##---------------------------------------------------------------------------##
## Publish all functions to the current shell.
##---------------------------------------------------------------------------##

#export -f whichall cleanemacs ssh1 scp1 npwd fixperms mpush mpull
#export -f loadandlistmodules findsymbol xe pkgdepends cuo

