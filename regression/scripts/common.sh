#!/bin/bash -l
## -*- Mode: sh -*-
##---------------------------------------------------------------------------##
## File  : regression/sripts/common.sh
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
## die           - exit with a message
## run           - echo a command and then run it.
## fn_exists     - return true if named bash function is defined
## establish_permissions - Change group to othello, dacodes or draco and change
##                 permissions to g+rwX,o-rwX
## machineName   - return a string to represent the current machine.
## osName        - return a string to represent the current machine's OS.
## flavor        - build a string that looks like fire-openmpi-2.0.2-intel-17.0.1
## selectscratchdir - find a scratch drive
## lookupppn     - return PE's per node.
## npes_build    - return PE's to be used for compiling.
## npes_test     - return PE's to be used for testing.
## install_verions - helper for doing releases (see release_toss2.sh)
## publish_release - helper for doing releases (see release_toss2.sh)
## allow_file_to_age - pause a program until a file is 'old'

##---------------------------------------------------------------------------##
## Helpful functions
##---------------------------------------------------------------------------##

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

function establish_permissions
{
  # Permissions - new files should be marked u+rwx,g+rwx,o+rx
  # Group is set to $1 or draco
  umask 0002
  if [[ `groups | grep -c othello` = 1 ]]; then
    install_group="othello"
    install_permissions="g+rwX,o-rwX"
  elif [[ `groups | grep -c dacodes` = 1 ]]; then
    install_group="dacodes"
    install_permissions="g+rwX,o-rwX"
  else
    install_group="draco"
    install_permissions="g+rwX,o-rwX"
  fi
  build_group="$USER"
  build_permissions="g+rwX,o-rwX"
}

# Logic taken from /usr/projects/hpcsoft/templates/header
function machineName
{
  sysName=${sysName="unknown"}
  if [[ -f /usr/projects/hpcsoft/utilities/bin/sys_name ]]; then
    sysName=`/usr/projects/hpcsoft/utilities/bin/sys_name`
  elif [[ -d /projects/darwin ]]; then
    sysName=darwin
  elif test -d /usr/gapps/jayenne; then
    sysName=sq
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
  elif [[ -d /projects/darwin ]]; then
    osName=darwin
  elif test -d /usr/gapps/jayenne; then
    osName=`uname -p`
  fi
  if [[ "$osName" == "unknown" ]]; then
    die "Unable to determine system OS, please edit scripts/common.sh."
  fi
  echo $osName
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
      OLDIFS=$IFS
      IFS=:
      # select modules that look like compilers
      unset compilermodules
      for module in $loadedmodules; do
        case $module in
          PrgEnv*)
            # Ingore PrgEnv matches.
            ;;
          intel/* | pgi/* | cray/* | gnu/* )
            tmp=`echo $module | sed -e 's%/%-%'`
            compilermodules="$tmp $compilermodules"
            ;;
        esac
      done
      IFS=$OLDIFS
      # pick the first compiler in the list
      compilerflavor=`echo $compilermodules | sed -e 's/ *//'`
      # append target if KNL
      if [[ `echo $CRAY_CPU_TARGET | grep -c knl` == 1 ]]; then
        compilerflavor+='-knl'
      fi
      ;;
    darwin*)
      if [[ $MPIARCH ]]; then
        if [[ $MPI_ROOT ]]; then
          LMPIVER=`echo $MPI_ROOT | sed -r 's%.*/([0-9]+)[.]([0-9]+)[.]([0-9]+).*%\1.\2.\3%'`
        else
          LMPIVER=''
        fi
        mpiflavor=$MPIARCH-$LMPIVER
      else
        mpiflavor="unknown"
      fi
      if [[ $LCOMPILER ]]; then
        compilerflavor=$LCOMPILER-$LCOMPILERVER
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
function selectscratchdir
{
  # if df is too old this command won't work correctly, use an alternate form.
  local scratchdirs=`df --output=pcent,target 2>&1 | grep -c unrecognized`
  if [[ $scratchdirs == 0 ]]; then
    scratchdirs=`df --output=pcent,target | grep scratch | grep -v netscratch | sort -g`
  else
    scratchdirs=`df -a 2> /dev/null | grep net/scratch | awk '{ print $4 " "$5 }' | sort -g`
    if ! [[ $scratchdirs ]]; then
      scratchdirs=`df -a 2> /dev/null | grep lustre/scratch | awk '{ print $4 " "$5 }' | sort -g`

    fi
  fi
  local odd=1
  for item in $scratchdirs; do
    # odd numbered items are disk's 'percent full'. They are ordered from least
    # used to most used.  Skip these values.
    if [[ $odd == 1 ]]; then
      odd=0
      continue
    else
      odd=1
    fi
    # if this location is good (must be able to write to this location), return
    # the path.
    mkdir -p $item/$USER &> /dev/null
    touch $item/$USER/selectscratchdir &> /dev/null
    if [[ -f $item/$USER/selectscratchdir ]]; then
      rm $item/$USER/selectscratchdir
      echo "$item"
      return
    fi
    # might need another directory level 'yellow'
    mkdir -p $item/yellow/$USER &> /dev/null
    touch $item/yellow/$USER/selectscratchdir &> /dev/null
    if [[ -f $item/yellow/$USER/selectscratchdir ]]; then
      rm $item/yellow/$USER/selectscratchdir
      echo "$item/yellow"
      return
    fi
  done

  # if no writable scratch directory is located, then also try netscratch;
  item=/netscratch/$USER
  mkdir -p $item &> /dev/null
  touch $item/selectscratchdir &> /dev/null
  if [[ -f $item/selectscratchdir ]]; then
    rm $item/selectscratchdir
    echo "$item"
    return
  fi

}

#------------------------------------------------------------------------------#
function lookupppn()
{
  # https://hpc.lanl.gov/index.php?q=summary_table
  local target="`uname -n | sed -e s/[.].*//`"
  local ppn=1
  case ${target} in
    ml* | pi* | wf* | lu* ) ppn=16 ;;
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
    fi* | ic* | sn* ) ppn=36 ;;
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
  elif [[ -f /proc/cpuinfo ]]; then
    # lscpu=`lscpu | grep "CPU(s):" | head -n 1 | awk '{ print $2 }'`
    np=`cat /proc/cpuinfo | grep -c processor`
  fi
  echo $np
}

function npes_test
{
  local np=1
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
  if test -z "${version}"; then
    echo "You must provide variable version."
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
  echo "# Begin release: $buildflavor/$version"
  echo "# ------------"
  echo

  # Create install directory
  install_dir="$install_prefix/$version"
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
  build_dir="$build_prefix/$version/${package:0:1}"

  # Purge any existing files before running cmake to configure the build directory.
  if test $config_step == 1; then
    if test -d ${build_dir}; then
      run "rm -rf ${build_dir}"
    fi
    run "mkdir -p $build_dir" || die "Could not create directory $build_dir."
  fi

  run "cd $build_dir"
  if test $config_step == 1; then
    run "cmake -DCMAKE_INSTALL_PREFIX=$install_dir \
             $options $CONFIG_EXTRA $source_dir" \
      || die "Could not configure in $build_dir from source at $source_dir"
  fi
  if test $build_step == 1; then
    run "make -j $build_pe -l $build_pe install"  \
      || die "Could not build code/tests in $build_dir"
  fi
  if test $test_step == 1; then
    case $version in
      *_nr)
        # only run tests that are safe in non-reproducible mode.
        run "ctest -L nr -j $test_pe" ;;
      *)
        # run all tests
        run "ctest -j $test_pe" ;;
    esac
  fi
  if ! test ${build_permissions:-notset} = "notset"; then
    run "chmod -R $build_permissions $build_dir"
  fi
}

##----------------------------------------------------------------------------##
## If $jobids is set, wait for those jobs to finish before setting
## groups and permissions.
function publish_release()
{
  echo " "
  echo "Waiting batch jobs to finish ..."
  echo "   Running jobs = $jobids"

  establish_permissions

  case `osName` in
    toss* | cle* ) SHOWQ=squeue ;;
    darwin| ppc64) SHOWQ=squeue ;;
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
function job_launch_sanity_checks()
{
  if [[ ! ${regdir} ]]; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'regdir' in the environment!"
    echo "printenv -> "
    printenv
    exit 1
  fi
  if [[ ! ${rscriptdir} ]]; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'rscriptdir' in the environment!"
    echo "printenv -> "
    printenv
    exit 1
  fi
  if [[ ! ${subproj} ]]; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'subproj' in the environment!"
    echo "printenv -> "
    printenv
    exit 1
  fi
  if [[ ! ${build_type} ]]; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'build_type' in the environment!"
    echo "printenv -> "
    printenv
    exit 1
  fi
  if [[ ! ${logdir} ]]; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'logdir' in the environment!"
    echo "printenv -> "
    printenv
    exit 1
  fi
  if [[ ! ${featurebranch} ]]; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'featurebranch' in the environment!"
    echo "printenv -> "
    printenv
  fi
}

#------------------------------------------------------------------------------#
# Job Launch Banner
function print_job_launch_banner()
{
echo "==========================================================================="
echo "$machine_name_long regression job launcher for ${subproj} - ${build_type} flavor."
echo "==========================================================================="
echo " "
echo "Environment:"
echo "   build_autodoc  = ${build_autodoc}"
echo "   build_type     = ${build_type}"
echo "   dashboard_type = ${dashboard_type}"
echo "   epdash         = $epdash"
if [[ ! ${extra_params} ]]; then
  echo "   extra_params   = none"
else
  echo "   extra_params   = ${extra_params}"
fi
if [[ ${featurebranch} ]]; then
  echo "   featurebranch  = ${featurebranch}"
fi
echo "   logdir         = ${logdir}"
echo "   logfile        = ${logfile}"
echo "   machine_name_long = $machine_name_long"
echo "   prdash         = $prdash"
echo "   projects       = \"${projects}\""
echo "   regdir         = ${regdir}"
echo "   regress_mode   = ${regress_mode}"
echo "   rscriptdir     = ${rscriptdir}"
echo "   scratchdir     = ${scratchdir}"
echo "   subproj        = ${subproj}"
echo " "
echo "   ${subproj}: dep_jobids = ${dep_jobids}"
echo " "
}

##----------------------------------------------------------------------------##
export die
export run
export establish_permissions
export machineName
export osName
export flavor
export selectscratchdir
export npes_build
export npes_test
export install_versions
export job_launch_sanity_checks

##----------------------------------------------------------------------------##
## End common.sh
##----------------------------------------------------------------------------##
