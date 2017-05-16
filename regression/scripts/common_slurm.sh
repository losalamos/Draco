#!/bin/bash -l

##---------------------------------------------------------------------------##
## Helpful functions
##---------------------------------------------------------------------------##

function die () { echo "ERROR: $1"; exit 1;}

function run () {
  echo "==> $1"
  if test ${dry_run:-no} = "no"; then eval $1; fi
}

fn_exists()
{
  type $1 2>/dev/null | grep -c 'is a function'
}

#----------------------------------------------------------------------#
function establish_permissions
{
  # Permissions - new files should be marked u+rwx,g+rwx,o+rx
  # Group is set to $1 or draco
  umask 0002
  if test `groups | grep -c othello` = 1; then
    install_group="othello"
    install_permissions="g+rwX,o-rwX"
  elif test `groups | grep -c dacodes` = 1; then
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
  if test -f /usr/projects/hpcsoft/utilities/bin/sys_name; then
    sysName=`/usr/projects/hpcsoft/utilities/bin/sys_name`
  elif test -d /projects/darwin; then
    sysName=darwin
  elif test -d /usr/gapps/jayenne; then
    sysName=sq
  fi
  if test "$sysName" = "unknown"; then
    echo "Unable to determine machine name, please edit scripts/common.sh."
    exit 1
  fi
  echo $sysName
}

# Logic taken from /usr/projects/hpcsoft/templates/header
function osName
{
  osName=${osName="unknown"}
  if test -f /usr/projects/hpcsoft/utilities/bin/sys_os; then
    osName=`/usr/projects/hpcsoft/utilities/bin/sys_os`
  elif test -d /projects/darwin; then
    osName=darwin
  elif test -d /usr/gapps/jayenne; then
    osName=`uname -p`
  fi
  if test "$osName" = "unknown"; then
    echo "Unable to determine system OS, please edit scripts/common.sh."
    exit 1
  fi
  echo $osName
}

function flavor
{
  platform=`machineName`
  os=`osName`
  case $os in
    toss*)
      if test -z $LMPI; then
        mpiflavor="unknown"
      else
        mpiflavor=$LMPI-$LMPIVER
      fi
      if test -z $LCOMPILER; then
        compilerflavor="unknown"
      else
        compilerflavor=$LCOMPILER-$LCOMPILERVER
      fi
      ;;
    cle*)
      if test -z $CRAY_MPICH2_VER; then
        mpiflavor="unknown"
      else
        mpiflavor=mpt-$CRAY_MPICH2_VER
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
      if test `echo $CRAY_CPU_TARGET | grep -c knl` == 1; then
        compilerflavor+='-knl'
      fi
      ;;
    darwin*)
      if test -z $MPIARCH; then
        mpiflavor="unknown"
      else
        if test -z $MPI_ROOT; then
          LMPIVER=''
        else
          LMPIVER=`echo $MPI_ROOT | sed -r 's%.*/([0-9]+)[.]([0-9]+)[.]([0-9]+).*%\1.\2.\3%'`
        fi
        mpiflavor=$MPIARCH-$LMPIVER
      fi
      if test -z $LCOMPILER; then
        compilerflavor="unknown"
      else
        compilerflavor=$LCOMPILER-$LCOMPILERVER
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
  esac
  echo $platform-$mpiflavor-$compilerflavor
}

function selectscratchdir
{
  # TOSS, CLE, BGQ, Darwin:
  toss2_yellow_scratchdirs="lustre/scratch1/yellow lustre/scratch2/yellow lustre/scratch3/yellow"
  toss2_red_scratchdirs="lustre/scratch3 lustre/scratch4"
  cray_yellow_scratchdirs="lustre/ttscratch1"
  cray_red_scratchdirs="lustre/trscratch1 lustre/trscratch2"
  bgq_scratchdirs="nfs/tmp2"
  scratchdirs="$toss2_yellow_scratchdirs $toss2_red_scratchdirs \
$cray_yellow_scratchdirs $cray_red_scratchdirs $bgq_scratchdirs \
usr/projects/draco/devs/releases"
  for dir in $scratchdirs; do
    mkdir -p /$dir/$USER &> /dev/null
    if test -x /$dir/$USER; then
      echo "$dir"
      return
    fi
  done
}

#------------------------------------------------------------------------------#
function lookupppn()
{
  # https://hpc.lanl.gov/index.php?q=summary_table
  local target="`uname -n | sed -e s/[.].*//`"
  local ppn=1
  case ${target} in
    ml* | pi* | wf* | lu* ) ppn=16 ;;
    tr-fe* | tr-login*) ppn=32 ;;
    tt-fe* | tt-login*)
      if [[ $CRAY_CPU_TARGET ]]; then
        if [[ $CRAY_CPU_TARGET == 'haswell' ]]; then
          ppn=32
        elif [[ $CRAY_CPU_TARGET == 'knl' ]]; then
          ppn=68
        fi
      else
        echo "ERROR: Expected CRAY_CPU_TARGET to be set in the environment."
        exit 1
      fi
      ;;
    fi* | ic* | sn* ) ppn=36 ;;
  esac
  echo $ppn
}

function npes_build
{
  local np=1
  slurm_num_cpus=`squeue -h -j $SLURM_JOB_ID -o "%C"`
  if [[ $slurm_num_cpus ]]; then
    np=$slurm_num_cpus
  elif ! test "${PBS_NP:-notset}" = "notset"; then
    np=${PBS_NP}
  elif ! test "${SLURM_NPROCS:-notset}" = "notset"; then
    np=${SLURM_NPROCS}
  elif ! test "${SLURM_CPUS_ON_NODE:-notset}" = "notset"; then
    np=${SLURM_CPUS_ON_NODE}
  elif ! test "${SLURM_TASKS_PER_NODE:-notset}" = "notset"; then
    np=${SLURM_CPUS_ON_NODE}
  elif test -f /proc/cpuinfo; then
    np=`cat /proc/cpuinfo | grep -c processor`
  fi
  echo $np
}

function npes_test
{
  local np=1
  slurm_num_cpus=`squeue -h -j $SLURM_JOB_ID -o "%C"`
  if [[ $slurm_num_cpus ]]; then
    np=$slurm_num_cpus
  elif ! test "${PBS_NP:-notset}" = "notset"; then
    np=${PBS_NP}
  elif ! test "${SLURM_NPROCS:-notset}" = "notset"; then
    np=${SLURM_NPROCS}
  elif ! test "${SLURM_CPUS_ON_NODE:-notset}" = "notset"; then
    np=${SLURM_CPUS_ON_NODE}
  elif ! test "${SLURM_TASKS_PER_NODE:-notset}" = "notset"; then
    np=${SLURM_CPUS_ON_NODE}
  elif test `uname -p` = "ppc"; then
    # sinfo --long --partition=pdebug (show limits)
    np=64
  elif test -f /proc/cpuinfo; then
    np=`cat /proc/cpuinfo | grep -c processor`
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
  if test -z ${build_pe}; then
    build_pe=`npes_build`
  fi
  if test -z ${test_pe}; then
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

    # dirs="$script_dir $source_prefix/source $source_prefix/logs"
    # for dir in $dirs; do
    #   if test -d $dir; then
    #     run "chgrp -R draco $dir"
    #     run "chmod $build_permissions $dir"
    #   fi
    # done

  fi
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
##----------------------------------------------------------------------------##
