#!/bin/bash -l
## -*- Mode: sh -*-
##---------------------------------------------------------------------------##
## File  : regression/sripts/common.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2019, Triad National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##
##
## Summary: Misc bash functions useful during development of code.

shopt -s extglob # Allow variable as case condition

##---------------------------------------------------------------------------##
## Helpful functions
##---------------------------------------------------------------------------##

function dracohelp()
{
echo -e "Bash functions defined by Draco:\n\n"
echo -e "Also try 'slurmhelp'\n"
echo "allow_file_to_age - pause a program until a file is 'old'"
echo "canonicalize_filename - standardize paths to files"
echo "cleanemacs    - remove ~ and .elc files."
echo "die           - exit with a message"
echo "dracoenv/rmdracoenv - load/unload the draco environment"
echo "establish_permissions - Change group to othello, dacodes or draco and change permissions to g+rwX,o-rwX"
echo "flavor        - build a string that looks like fire-openmpi-2.0.2-intel-17.0.1"
echo "fn_exists     - return true if named bash function is defined"
echo "install_verions - helper for doing releases (see release_toss2.sh)"
echo "lookupppn     - return PE's per node."
echo "machineName   - return a string to represent the current machine."
echo "matches_extension - Does the provided filename have a matching extension?"
echo "npes_build    - return PE's to be used for compiling."
echo "npes_test     - return PE's to be used for testing."
echo "osName        - return a string to represent the current machine's OS."
echo "proxy         - toggle the status of the LANL proxy variables."
echo "publish_release - helper for doing releases (see release_toss2.sh)"
echo "qrm           - quick rm for directories located in lustre scratch spaces."
echo "rdde          - more agressive reset of the draco environment."
echo "run           - echo a command and then run it."
echo "selectscratchdir - find a scratch drive"
echo "version_gt    - compare versions"
echo "whichall      - Find all matchs."
echo "xfstatus      - print status 'transfered' files."
echo -e "\nUse 'type <function>' to print the full content of any function.\n"
}

#------------------------------------------------------------------------------#
# Function Definitions
#------------------------------------------------------------------------------#

# Print an error message and exit.
# e.g.: cd $dir || die "can't change dir to $dir".
function die () { echo " "; echo "FATAL ERROR: $1"; exit 1;}

# Echo a command and then run it.
function run ()
{
  echo "==> $1"; if test ${dry_run:-no} = "no"; then eval $1; fi
}

# Return 0 if provided name is a bash function.
function fn_exists ()
{
  type $1 2>/dev/null | grep -c 'is a function'
}

#----------------------------------------------------------------------#
# The script starts here
#----------------------------------------------------------------------#

# Logic taken from /usr/projects/hpcsoft/templates/header
function machineName
{
  sysName=${sysName="unknown"}
  if [[ -f /usr/projects/hpcsoft/utilities/bin/sys_name ]]; then
    sysName=`/usr/projects/hpcsoft/utilities/bin/sys_name`
  elif [[ -d /projects/darwin ]]; then
    sysName=darwin
  elif [[ -d /usr/gapps/jayenne ]]; then
    case `uname -n` in
      sq* | sequoia*) sysName=sq ;;
      rzansel*) sysName=rzansel ;;
      rzmanta*) sysName=rzmanta ;;
      sierra*) sysName=sierra ;;
    esac
  fi
  if [[ "$sysName" == "unknown" ]]; then
    die "Unable to determine machine name, please edit scripts/common.sh."
  fi
  echo $sysName
}

# Logic taken from /usr/projects/hpcsoft/templates/header
function osName
{
  osName=${osName="unknown"}
  if [[ -f /usr/projects/hpcsoft/utilities/bin/sys_os ]]; then
    osName=`/usr/projects/hpcsoft/utilities/bin/sys_os`
    if [[ $? != 0 ]]; then
      osName="unknown"
    fi
  fi
  if [[ ${osName} == "unknown" ]]; then
    if [[ -d /projects/darwin ]] ; then
      osName=darwin
      if [[ -f /projects/draco/vendors/bin/target_arch ]]; then
        osName=darwin-`/projects/draco/vendors/bin/target_arch`
      fi
    elif [[ -d /usr/gapps/jayenne ]]; then
      osName=`uname -p`
    fi
  fi
  if [[ "$osName" == "unknown" ]]; then
    die "Unable to determine system OS, please edit scripts/common.sh."
  fi
  echo $osName
}

function machineFamily
{
  machfam=`uname -p`
  case `osName` in
    darwin*) machfam=`osName` ;;
    cle*)
      case `uname -n | sed -e 's/[.].*//'` in
        tt* | tr* ) machfam="ats1" ;;
        cp* | th* ) machfam="cray" ;;
      esac
      ;;
    toss3) machfam="cts1" ;;
    ppc64le)
     case `uname -n | sed -e 's/[.].*//'` in
       rzans* | sierra* ) machfam="ats2" ;;
     esac
     ;;
  esac
  echo $machfam
}

#------------------------------------------------------------------------------#
# Generates a string of the form <platform>-<mpi+ver>-<compiler+ver>
function flavor
{
  platform=`machineName`
  os=`osName`
  case $os in
    toss*)
      if [[ $LMPI ]]; then
        mpiflavor=$LMPI-$LMPIVER
      else
        mpiflavor="unknown"
      fi
      if [[ $LCOMPILER ]]; then
        compilerflavor=$LCOMPILER-$LCOMPILERVER
      else
        compilerflavor="unknown"
      fi
      ;;
    cle*)
      if [[ $CRAY_MPICH2_VER ]]; then
        mpiflavor=mpt-$CRAY_MPICH2_VER
      else
        mpiflavor="unknown"
      fi
      # Try to determine the loaded compiler
      loadedmodules=`echo $LOADEDMODULES`
      # select modules that look like compilers
      unset compilermodules
      for module in ${loadedmodules//:/ }; do
        case $module in
          PrgEnv*)
            # Ingore PrgEnv matches.
            ;;
          cce/* | cray/* | gcc/* | gnu/* | intel/* | pgi/* )
            tmp=`echo $module | sed -e 's%/%-%'`
            compilermodules="$tmp $compilermodules"
            ;;
        esac
      done
      # pick the first compiler in the list
      compilerflavor=`echo $compilermodules | sed -e 's/ *//'`
      # append target if KNL
      if [[ `echo $CRAY_CPU_TARGET | grep -c knl` == 1 ]]; then
        compilerflavor+='-knl'
      fi
      ;;
    darwin*)
      platform=$os # e.g.: darwin-power9
      if [[ $MPIARCH ]]; then
        if [[ $MPIRUN ]]; then
          LMPIVER=`echo $MPIRUN | sed -r 's%.*/([0-9]+)[.]([0-9]+)[.]([0-9]+).*%\1.\2.\3%'`
        else
          LMPIVER=''
        fi
        mpiflavor=$MPIARCH-$LMPIVER
      else
        mpiflavor="unknown"
      fi
      if [[ $CC ]]; then
        if [[ $CC =~ "gcc" ]]; then
          compilerflavor=gnu-`$CC --version | head -n 1 | sed -r 's%.* %%'`
        fi
      else
        compilerflavor="unknown"
      fi
      ;;
    ppc64)
      # more /bgsys/drivers/V1R2M3/ppc64/comm/include/mpi.h
      # | grep MPI_VERSION    ==> 2 ==> (mpich2)
      # | grep MPICH2_VERSION ==> 1.5
      mpiflavor="mpich2-1.5"

      case $CC in
      *gcc*)
          LCOMPILER=gnu
          LCOMPILERVER=`$CC --version | head -n 1 | sed -e 's/.*\([0-9][.][0-9][.][0-9]\)/\1/'`
          compilerflavor=$LCOMPILER-$LCOMPILERVER
          ;;
      *xlc*)
          LCOMPILER=ibm
          LCOMPILERVER=`$CC -V | head -n 1 | sed -e 's/.*[/]\([0-9]\+\.[0-9]\).*/\1/'`
          compilerflavor=$LCOMPILER-$LCOMPILERVER
          ;;
      *)
          compiler_flavor=unknown-unknown ;;
      esac
      ;;

    ppc64le)
      if [[ -n $OPAL_PREFIX ]]; then
        mpiflavor=`echo $OPAL_PREFIX | sed -e 's%.*/%%'`
      elif [[ -n $LMOD_MPI_NAME ]]; then
        mpiver=`echo $LMOD_MPI_VERSION | sed -e 's%-.*%%'`
        mpiflavor=$LMOD_MPI_NAME-$mpiver
      else
        mpiflavor=unknown
      fi
      case $CC in
      *gcc*)
          LCOMPILER=gnu
          LCOMPILERVER=`$CC --version | head -n 1 | sed -e 's/.*\([0-9][.][0-9][.][0-9]\).*/\1/'`
          compilerflavor=$LCOMPILER-$LCOMPILERVER
          ;;
      *xlc*)
          LCOMPILER=ibm
          LCOMPILERVER=`$CC --version | head -n 1 | sed -e 's/.*V\([0-9][0-9][.][0-9][.][0-9]\) .*/\1/'`
          compilerflavor=$LCOMPILER-$LCOMPILERVER
          ;;
      *) compiler_flavor=unknown-unknown ;;
      esac
      ;;

    *)
      # CCS-NET machines or generic Linux?
      if [[ $MPI_NAME ]]; then
        mpiflavor=$MPI_NAME-$MPI_VERSION
      else
        mpiflavor="unknown"
      fi
      if [[ $LCOMPILER ]]; then
        compilerflavor=$LCOMPILER-$LCOMPILERVER
      else
        compilerflavor="unknown"
      fi
      ;;
  esac
  echo $platform-$mpiflavor-$compilerflavor
}

#------------------------------------------------------------------------------#
# returns a path to a directory
function selectscratchdir ()
{
    local scratchdirs=`df --output=pcent,target 2>&1 | grep -c unrecognized`;
    if [[ $scratchdirs == 0 ]]; then
        scratchdirs=`df --output=pcent,target 2> /dev/null | grep scratch | grep -v netscratch | sort -g`;
        if [[ -z "$scratchdirs" ]]; then
            scratchdirs=`df --output=pcent,target 2> /dev/null | grep workspace | sort -g`;
        fi;
    else
        scratchdirs=`df -a 2> /dev/null | grep net/scratch | awk '{ print $4 " "$5 }' | sort -g`;
        if [[ -z "$scratchdirs" ]]; then
            scratchdirs=`df -a 2> /dev/null | grep lustre/scratch | awk '{ print $4 " "$5 }' | sort -g`;
        fi;
    fi;
    for item in $scratchdirs;
    do
        if [[ `echo $item | grep -c %` == 1 ]]; then
            continue;
        fi;
        if [[ -d $item/users ]]; then
            item2="${item}/users";
            mkdir -p $item2/$USER &>/dev/null;
            if [[ -w $item2/$USER ]]; then
                echo "$item2";
                return;
            fi;
        fi;
        mkdir -p $item/$USER &>/dev/null;
        if [[ -w $item/$USER ]]; then
            echo "$item";
            return;
        fi;
        mkdir -p $item/yellow/$USER &>/dev/null;
        if [[ -w $item/yellow/$USER ]]; then
            echo "$item/yellow";
            return;
        fi;
    done;
    item=/netscratch/$USER;
    mkdir -p $item &>/dev/null;
    if [[ -w $item ]]; then
        echo "$item";
        return;
    fi
}

#------------------------------------------------------------------------------#
function lookupppn()
{
  # https://hpc.lanl.gov/index.php?q=summary_table
  local target="`uname -n | sed -e s/[.].*//`"
  local ppn=1
  case ${target} in
    t[rt]-fe* | t[rt]-login*)
      if [[ $CRAY_CPU_TARGET == "haswell" ]]; then
          ppn=32
      elif [[ $CRAY_CPU_TARGET == "knl" ]]; then
        ppn=68
      else
        echo "ERROR: Expected CRAY_CPU_TARGET to be set in the environment."
        exit 1
      fi
      ;;
    fi* | ic* | sn* | cy* )         ppn=36 ;;
    rzansel* | rzmanta* | sierra* ) ppn=40 ;;
    *) ppn=`cat /proc/cpuinfo | grep -c processor` ;;
  esac
  echo $ppn
}

function npes_build
{
  local np=1
  if [[ ${PBS_NP} ]]; then
    np=${PBS_NP}
  elif [[ ${SLURM_NPROCS} ]]; then
    np=${SLURM_NPROCS}
  elif [[  ${SLURM_CPUS_ON_NODE} ]]; then
    np=${SLURM_CPUS_ON_NODE}
  elif [[ ${SLURM_TASKS_PER_NODE} ]]; then
    np=${SLURM_TSKS_PER_NODE}
  elif [[ ${LSB_DJOB_NUMPROC} ]]; then
    np=${LSB_DJOB_NUMPROC}
  elif [[ -f /proc/cpuinfo ]]; then
    # lscpu=`lscpu | grep "CPU(s):" | head -n 1 | awk '{ print $2 }'`
    np=`cat /proc/cpuinfo | grep -c processor`
  fi
  echo $np
}

function npes_test
{
  # assume no parallel testing capability.
  local np=1

  # allow specialization per machine (if needed)
  local target="`uname -n | sed -e s/[.].*//`"
  case ${target} in

    # current LSF only allows one executable to run under 'jsrun' at a time.
    # rzansel* | rzmanta* | sierra*) ppn=1 ;;

    *)

      # use lscpu if it is available.
      if ! [[ `which lscpu 2>/dev/null` == 0 ]]; then
        # number of cores per socket
        local cps=`lscpu | grep "^Core(s)" | awk '{ print $4 }'`
        # number of sockets
        local ns=`lscpu | grep "^Socket(s):" | awk '{ print $2 }'`
        np=`expr $cps \* $ns`

      else

        if [[ ${PBS_NP} ]]; then
          np=${PBS_NP}
        elif [[ ${SLURM_NPROCS} ]]; then
          np=${SLURM_NPROCS}
        elif [[  ${SLURM_CPUS_ON_NODE} ]]; then
          np=${SLURM_CPUS_ON_NODE}
        elif [[ ${SLURM_TASKS_PER_NODE} ]]; then
          np=${SLURM_TSKS_PER_NODE}
        elif [[ `uname -p` == "ppc" ]]; then
          # sinfo --long --partition=pdebug (show limits)
          np=64
        elif [[ -f /proc/cpuinfo ]]; then
          # lscpu=`lscpu | grep "CPU(s):" | head -n 1 | awk '{ print $2 }'`
          np=`cat /proc/cpuinfo | grep -c processor`
        fi

      fi
      ;;
  esac
  echo $np
}

##---------------------------------------------------------------------------##
## Configure, Build and Run the Tests
##---------------------------------------------------------------------------##

function install_versions
{
  local config_step=0
  local build_step=0
  local test_step=0

  if test -z "${steps}"; then
    echo "You must provide variable steps."
    echo "E.g.: steps=\"configure build test\""
    return
  else
    for s in $steps; do
      case $s in
        config) config_step=1 ;;
        build)  build_step=1  ;;
        test)   test_step=1   ;;
      esac
    done
  fi
  if test -z ${buildflavor}; then
    buildflavor=`flavor`
  fi
  if test -z "${rttversion}"; then
    echo "You must provide variable rttversion."
    # echo "E.g.: VERSIONS=( \"debug\" \"opt\" )"
    return
  fi
  if test -z "${options}"; then
    echo "You must provide variable option."
    #echo "E.g.: OPTIONS=("
    #echo "\"$OPTIMIZE_OFF $LOGGING_OFF\""
    #echo "\"$OPTIMIZE_ON  $LOGGING_OFF\""
    #echo ")"
    return
  fi
  if test -z "${package}"; then
    echo "You must provide variable package."
    echo "E.g.: package=\"draco\""
    return
  fi
  if test -z "${install_prefix}"; then
    echo "You must provide variable install_prefix."
    echo "E.g.: install_prefix=/usr/projects/draco/$pdir/$buildflavor"
    return
  fi
  if ! [[ ${build_pe} ]]; then
    build_pe=`npes_build`
  fi
  if ! [[ ${test_pe} ]]; then
    test_pe=`npes_test`
  fi

  # Echo environment before we start:
  echo
  echo
  echo "# Environment"
  echo "# ------------"
  echo
  run "module list"
  run "printenv"
  echo "---"
  echo "Environment size = `printenv | wc -c`"

  echo
  echo
  echo "# Begin release: $buildflavor/$rttversion"
  echo "# ------------"
  echo

  # Create install directory
  install_dir="$install_prefix/$rttversion"
  if ! test -d $install_dir; then
    run "mkdir -p $install_dir" || die "Could not create $install_dir"
  fi

  # try to locate the souce directory
  if test -f $source_prefix/source/ChangeLog; then
    source_dir=$source_prefix/source
  else
    local possible_source_dirs=`/bin/ls -1 $source_prefix/source`
    for dir in "$possible_source_dirs"; do
      if test -f $source_prefix/source/$dir/ChangeLog; then
        source_dir=$source_prefix/source/$dir
        break
      fi
    done
  fi
  if ! test -f $source_dir/CMakeLists.txt; then
    echo "Could not find sources. Tried looking at $source_prefix/source/"
    exit 1
  fi
  # source_dir="$source_prefix/source/$package"
  build_dir="$build_prefix/$rttversion/${package:0:1}"

  # Purge any existing files before running cmake to configure the build
  # directory.
  if test $config_step == 1; then
    if test -d ${build_dir}; then
      run "rm -rf ${build_dir}"
    fi
    run "mkdir -p $build_dir" || die "Could not create directory $build_dir."
  fi

  run "cd $build_dir"
  if test $config_step == 1; then
    run "cmake -DCMAKE_INSTALL_PREFIX=$install_dir $options $CONFIG_EXTRA $source_dir" \
      || die "Could not configure in $build_dir from source at $source_dir"
  fi
  if test $build_step == 1; then
    run "make -j $build_pe -l $build_pe install"  \
      || die "Could not build code/tests in $build_dir"
  fi
  if test $test_step == 1; then
    # run all tests
    run "ctest -j $test_pe --output-on-failure"
    if [[ $? != 0 ]]; then
      run "ctest -j $test_pe --output-on-failure --rerun-failed"
    fi
  fi
  wait
  if ! test ${build_permissions:-notset} = "notset"; then
    run "chmod -R $build_permissions $build_dir"
  fi
}

##----------------------------------------------------------------------------##
## If $jobids is set, wait for those jobs to finish before setting groups and
## permissions.
function publish_release()
{
  echo " "
  echo "Waiting batch jobs to finish ..."
  echo "   Running jobs = $jobids"

  case `osName` in
    toss* | cle* ) SHOWQ=squeue ;;
    darwin| ppc64) SHOWQ=squeue ;;
    ppc64le )      SHOWQ=bjobs  ;;
  esac

  # wait for jobs to finish
  for jobid in $jobids; do
    while test `${SHOWQ} | grep -c $jobid` -gt 0; do
      ${SHOWQ} | grep $jobid
      sleep 5m
    done
    echo "   Job $jobid is complete."
  done

  echo " "
  echo "Updating file permissions ..."
  echo " "

  if ! test -z $install_permissions; then
    # Set access to top level install dir.
    if test -d $install_prefix; then
      run "chgrp -R ${install_group} $source_prefix"
      run "chmod -R $install_permissions $source_prefix"
      run "find $source_prefix -type d -exec chmod g+s {} +"
    fi
  fi
}

#------------------------------------------------------------------------------#
# Pause until the 'last modified' timestamp of file $1 to be $2 seconds old.
function allow_file_to_age
{
  if [[ ! $2 ]]; then
    echo "ERROR: This function requires two arguments: a filename and an age value (sec)."
    exit 1
  fi

  # If file does not exist, no need to wait.
  if [[ ! -f $1 ]]; then
    return
  fi

  # assume file was last modified 0 seconds ago.
  local timediff=0

  # If no changes for $2 seconds, continue
  # else, wait until until file, $1, hasn't been touched for $2 seconds.
  local print_message=1
  while [[ $timediff -lt $2 ]]; do
    eval "$(date +'now=%s')"
    local pr_last_check=$(date +%s -r $1)
    local timediff=$(expr $now - $pr_last_check)
    local timeleft=$(expr $2 - $timediff)
    if [[ $timeleft -gt 0 ]]; then
      if [[ $print_message == 1 ]]; then
        echo "The log file $1 was recently modified."
        echo "To avoid colliding with another running test we are waiting"
        print_message=0
      fi
      echo "... $timeleft seconds"
    fi
    sleep 30s
  done
}

#------------------------------------------------------------------------------#
# Ensure environment is reasonable
# Called from <machine>-job-launch.sh
# function job_launch_sanity_checks()
# {
#   if [[ ! ${regdir} ]]; then
#     echo "FATAL ERROR in ${scriptname}: You did not set 'regdir' in the environment!"
#     echo "printenv -> "
#     printenv
#     exit 1
#   fi
#   if [[ ! ${rscriptdir} ]]; then
#     echo "FATAL ERROR in ${scriptname}: You did not set 'rscriptdir' in the environment!"
#     echo "printenv -> "
#     printenv
#     exit 1
#   fi
#   if [[ ! ${subproj} ]]; then
#     echo "FATAL ERROR in ${scriptname}: You did not set 'subproj' in the environment!"
#     echo "printenv -> "
#     printenv
#     exit 1
#   fi
#   if [[ ! ${build_type} ]]; then
#     echo "FATAL ERROR in ${scriptname}: You did not set 'build_type' in the environment!"
#     echo "printenv -> "
#     printenv
#     exit 1
#   fi
#   if [[ ! ${logdir} ]]; then
#     echo "FATAL ERROR in ${scriptname}: You did not set 'logdir' in the environment!"
#     echo "printenv -> "
#     printenv
#     exit 1
#   fi
#   if [[ ! ${featurebranch} ]]; then
#     echo "FATAL ERROR in ${scriptname}: You did not set 'featurebranch' in the environment!"
#     echo "printenv -> "
#     printenv
#   fi
# }

#------------------------------------------------------------------------------#
# Job Launch Banner
# function print_job_launch_banner()
# {
# echo "==========================================================================="
# echo "$machine_name_long regression job launcher for ${subproj} - ${build_type} flavor."
# echo "==========================================================================="
# echo " "
# echo "Environment:"
# echo "   build_autodoc  = ${build_autodoc}"
# echo "   build_type     = ${build_type}"
# echo "   dashboard_type = ${dashboard_type}"
# echo "   epdash         = $epdash"
# if [[ ! ${extra_params} ]]; then
#   echo "   extra_params   = none"
# else
#   echo "   extra_params   = ${extra_params}"
# fi
# if [[ ${featurebranch} ]]; then
#   echo "   featurebranch  = ${featurebranch}"
# fi
# echo "   logdir         = ${logdir}"
# echo "   logfile        = ${logfile}"
# echo "   machine_name_long = $machine_name_long"
# echo "   prdash         = $prdash"
# echo "   projects       = \"${projects}\""
# echo "   regdir         = ${regdir}"
# echo "   regress_mode   = ${regress_mode}"
# echo "   rscriptdir     = ${rscriptdir}"
# echo "   scratchdir     = ${scratchdir}"
# echo "   subproj        = ${subproj}"
# echo " "
# echo "   ${subproj}: dep_jobids = ${dep_jobids}"
# echo " "
# }

##----------------------------------------------------------------------------##
# Canonicalize_filename:
# See http://stackoverflow.com/questions/1055671/how-can-i-get-the-behavior-of-gnus-readlink-f-on-a-mac
##----------------------------------------------------------------------------##
canonicalize_filename () {
    local target_file=$1
    local physical_directory=""
    local result=""

    # Need to restore the working directory after work.
    pushd `pwd` > /dev/null

    cd "$(dirname "$target_file")"
    target_file=`basename $target_file`

    # Iterate down a (possible) chain of symlinks
    while [ -L "$target_file" ]
    do
        target_file=$(readlink "$target_file")
        cd "$(dirname "$target_file")"
        target_file=$(basename "$target_file")
    done

    # Compute the canonicalized name by finding the physical path
    # for the directory we're in and appending the target file.
    physical_directory=`pwd -P`
    result="$physical_directory"/"$target_file"

    # restore the working directory after work.
    popd > /dev/null

    echo "$result"
}

##----------------------------------------------------------------------------##
# Compare versions
#
# Example:
#
# EMACSVER=`emacs --version | head -n 1 | sed -e 's/.*Emacs //'`
# if `version_gt "24.0.0" $EMACSVER` ; then echo yes; fi
##----------------------------------------------------------------------------##
function version_gt()
{
  test "$(printf '%s\n' "$@" | sort -V | head -n 1)" != "$1";
}

##----------------------------------------------------------------------------##
# check whether the given file matches any of the set extensions
#
# Example:
#
# FILE_EXTS=".f90 .F90"
# FILE_ENDINGS="_f.h _f77.h _f90.h"
# for file in $modified_files; do
#   if ! matches_extension $file; then
#      continue;
#   fi
#   <do stuff>
# done
##----------------------------------------------------------------------------##
matches_extension() {
    local filename=$(basename "$1")
    local extension=".${filename##*.}"
    local end
    local ext
    for end in $FILE_ENDINGS; do [[ "$filename" == *"$end" ]] && return 1; done
    for ext in $FILE_EXTS; do [[ "$ext" == "$extension" ]] && return 0; done
    return 1
}

##----------------------------------------------------------------------------##
# relcreatesymlinks
#
# Creates symlinks to an installation to satisfy requests from clients.  For
# example:
#
#   cts1-openmpi-2.1.2-gcc-7.4.0    -> snow-openmpi-2.1.2-gcc-7.4.0
#   grizzly-openmpi-2.1.2-gcc-7.4.0 -> snow-openmpi-2.1.2-gcc-7.4.0
#
# This function is designed to be called from release-<machine>.sh
# The following variables should be set in the environment:
#
# - source_prefix :: e.g. /usr/projects/draco/draco-NN_NN_NN
# - evironments   :: a list of bash function names. These functiosn should
#                    establish build environments (see cts1-env.sh, etc.)
# - siblings      :: a space delimited list of machine names that are considered
#                    to be equivalent (e.g.: "fire ice cyclone")
#
##----------------------------------------------------------------------------##
function relcreatesymlinks() {

  OLDPWD=`pwd`
  machfam=`machineFamily`

  if ! [[ -d $source_prefix ]]; then
    die "relcreatesymlinks:: Expected variable source_prefix to be set"
  fi
  if ! [[ $siblings ]]; then
    die "relcreatesymlinks:: Expected variable siblings to be set"
  fi
  if ! [[ $environments ]]; then
    die "relcreatesymlinks:: Expected variable environments to be set"
  fi

  run "cd $source_prefix"

  for env in $environments; do

    if ! [[ `fn_exists $env` ]]; then
      die "relcreatesymlinks:: Attempted to load environment $env, but this function is not defined".
    fi
    # establish environment
    $env
    buildflavor=`flavor`

    # create symlinks
    me=`echo $buildflavor | sed -e 's/-.*//'`
    familymatch=`echo $siblings | sed -e 's/[ ]/|/g'`
    case $me in
      @($familymatch) )
        linkname=`echo $buildflavor | sed -e "s/${me}/${machfam}/"`
        run "ln -s $buildflavor $linkname"
        for m in $siblings; do
          linkname=`echo $buildflavor | sed -e "s/${me}/${m}/"`
          if ! [[ -d $source_prefix/$linkname ]]; then
            run "ln -s $buildflavor $linkname"
          fi
        done
        ;;
    esac
  done

  run "cd $OLDPWD"
}

##----------------------------------------------------------------------------##
export canonicalize_filename
export die
export flavor
export install_versions
export job_launch_sanity_checks
export machineName
export matches_extension
export npes_build
export npes_test
export osName
export relcreatesymlinks
export run
export selectscratchdir
export version_gt

##----------------------------------------------------------------------------##
## End common.sh
##----------------------------------------------------------------------------##
