#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/regression-master.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2019, Triad National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# Use:
# - Call from crontab using
#   <path>/regression-master.sh [options]

##---------------------------------------------------------------------------##
## Environment
##---------------------------------------------------------------------------##

# Because of this next 'exec sg' command, the crontab must escape double quotes
# to keep space delimited options together.  Something like:
# 00 06 * * 0-6 /scratch/regress/draco/regression/regression-master.sh -r -b Debug -d
# Nightly -p \"draco jayenne capsaicin core\" -e clang

# switch to group 'ccsrad' and set umask
if [[ $(id -gn) != ccsrad ]]; then
  exec sg ccsrad "$0 $*"
fi
umask 0007

# Enable job control
set -m

# Allow variable as case condition
shopt -s extglob

# load some common bash functions
export rscriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [[ -f $rscriptdir/scripts/common.sh ]]; then
  source $rscriptdir/scripts/common.sh
else
  echo " "
  echo "FATAL ERROR: Unable to locate Draco's bash functions: "
  echo "   looking for .../regression/scripts/common.sh"
  echo "   searched rscriptdir = $rscriptdir"
  exit 1
fi

# Host based variables
platform_extra_params=`echo $platform_extra_params | sed -e 's/ / | /g'`
export host=`uname -n | sed -e 's/[.].*//g'`
case $host in
  ba*|gr*|sn*) source $rscriptdir/cts1-options.sh ;;
  ccscs*)  source $rscriptdir/ccscs-options.sh ;;
  tt*)     source $rscriptdir/tt-options.sh ;;
  *)
    echo "FATAL ERROR: I don't know how to run regression on host = ${host}."
    print_use;  exit 1 ;;
esac

##---------------------------------------------------------------------------##
## Support functions
##---------------------------------------------------------------------------##

print_use()
{
  platform_extra_params=`echo $platform_extra_params | sed -e 's/ / | /g'`

  echo " "
  echo "Usage: ${0##*/} -b [Release|Debug] -d [Experimental|Nightly|Continuous]"
  echo "       -h -p [\"draco jayenne capsaicin core\"] -r"
  echo "       -f <git branch name> -a"
  echo "       -e <platform specific options>"
  echo " "
  echo "All arguments are optional,  The first value listed is the default value."
  echo "   -h    help           prints this message and exits."
  echo "   -r    regress        nightly regression mode."
  echo " "
  echo "   -a    build autodoc"
  echo "   -b    build-type     = { Debug, Release }"
  echo "   -d    dashboard type = { Experimental, Nightly, Continuous }"
  echo "   -f    git feature branch, default=\"develop develop\""
  echo "         common: 'develop pr42'"
  echo "         requires one string per project listed in option -p"
  echo "   -p    project names  = { draco, jayenne, capsaicin, core }"
  echo "                          This is a space delimited list within double quotes."
  echo "   -e    extra params   = { none | $platform_extra_params }"
  echo " "
  echo "Example:"
  echo "./regression-master.sh -b Release -d Nightly -p \"draco jayenne capsaicin, core\""
  echo " "
  echo "If no arguments are provided, this script will run"
  echo "   /regression-master.sh -b Debug -d Experimental -p \"draco\" -e none"
  echo " "
}

##---------------------------------------------------------------------------##
## Default values
##---------------------------------------------------------------------------##
build_autodoc="off"
build_type=Debug
dashboard_type=Experimental
projects="draco"
extra_params=""
regress_mode="off"
epdash=""
scratchdir=`selectscratchdir`

# Default to using GitHub for Draco and gitlab.lanl.gov for Jayenne
prdash="-"
unset nfb

##---------------------------------------------------------------------------##
## Command options
##---------------------------------------------------------------------------##

while getopts ":ab:d:e:f:ghp:r" opt; do
case $opt in
a)  build_autodoc="on";;
b)  build_type=$OPTARG ;;
d)  dashboard_type=$OPTARG ;;
e)  extra_params=$OPTARG
    epdash="-";;
f)  featurebranches=$OPTARG
    nfb=`echo $featurebranches | wc -w` ;;
h)  print_use; exit 0 ;;
p)  projects=$OPTARG ;;
r)  regress_mode="on" ;;
\?) echo "" ;echo "invalid option: -$OPTARG"; print_use; exit 1 ;;
:)  echo "" ;echo "option -$OPTARG requires an argument."; print_use; exit 1 ;;
esac
done
if [[ ${nfb} ]]; then
  # manually selecting feature branches -> must provide the same number of
  # feature branches as projects.
  if [[ ! `echo $projects | wc -w` == $nfb ]]; then
    echo "Error: You must provide the same number of feature branches as the number of"
    echo "projects specified. For example:"
    echo "    -p \"draco jayenne\" -f \"develop pr42\""
    exit 1
  fi
else
  # default: use 'develop' for all git branches.
  featurebranches=''
  for p in $projects; do
    featurebranches+="develop "
  done
fi

##---------------------------------------------------------------------------##
## Sanity Checks for input
##---------------------------------------------------------------------------##

case ${build_type} in
"Debug" | "Release" ) # known $build_type, continue
    ;;
*)  echo "" ;echo "FATAL ERROR: unsupported build_type (-b) = ${build_type}"
    print_use; exit 1 ;;
esac

case ${dashboard_type} in
Nightly | Experimental | Continuous) # known dashboard_type, continue
    ;;
*)  echo "" ;echo "FATAL ERROR: unknown dashboard_type (-d) = ${dashboard_type}"
    print_use; exit 1 ;;
esac

for proj in ${projects}; do
   case $proj in
   draco | jayenne | core | trt | npt ) # known projects, continue
      ;;
   *)  echo "" ;echo "FATAL ERROR: unknown project name (-p) = ${proj}"
       print_use; exit 1 ;;
   esac
done

case $regress_mode in
on)
    if ! [[ ${USER} == "kellyt" ]]; then
      echo "You are not authorized to use option '-r'."
      exit 1
    fi
    if [[ -d /usr/projects/jayenne/regress ]]; then
      regdir=/usr/projects/jayenne/regress
    else
      regdir=/scratch/regress
    fi
    logdir=$regdir/logs
    ;;
off)
    regdir="$scratchdir/$USER"
    logdir="$regdir/logs"
    ;;
*)  echo "" ;echo "FATAL ERROR: value of regress_mode=$regress_mode is incorrect."
    exit 1 ;;
esac

# Extra parameters valid for this machine?
if [[ ${extra_params} ]]; then
  for ep in $extra_params; do
    case $ep in
      none) extra_params=""; epdash="" ;;
      @($pem_match) )
      # known, continue
      ;;
      *)
        echo "" ;echo "FATAL ERROR: unknown extra params (-e) = ${extra_params}"
        print_use; exit 1
        ;;
    esac
  done
fi
# sort entries alphebetically to avoid knl_perfbench != perfbench_knl
extra_params_sort_safe=`echo $extra_params | xargs -n 1 | sort -u | xargs | sed -e 's/ /_/g'`

##---------------------------------------------------------------------------##
## Main
##---------------------------------------------------------------------------##

# Ensure log dir exists.
mkdir -p $logdir || die "Could not create a directory for log files."

# Redirect output to logfile.
timestamp=`date +%Y%m%d-%H%M`
logfile=$logdir/${machine_name_short}-${build_type}-master-$timestamp.log
case $regress_mode in
  off) echo "Redirecting output to $logfile" ;;
esac
exec > $logfile
exec 2>&1

##---------------------------------------------------------------------------##
## Export environment
##---------------------------------------------------------------------------##
export build_autodoc build_type dashboard_type epdash extra_params
export extra_params_sort_safe featurebranches logdir prdash machine_name_short
export regdir regress_mode scratchdir

##---------------------------------------------------------------------------##
# Banner
##---------------------------------------------------------------------------##

echo "==========================================================================="
echo "regression-master.sh: Regression for $machine_name_long ($machine_name_short)"
#echo "Build: ${build_type}     Extra Params: $extra_params"
date
echo "==========================================================================="
echo " "
echo "Host: $host"
echo " "
echo "Environment:"
echo "   build_autodoc  = $build_autodoc"
echo "   build_type     = ${build_type}"
echo "   dashboard_type = ${dashboard_type}"
echo "   epdash         = $epdash"
echo "   extra_params   = \"${extra_params}\" (${extra_params_sort_safe})"
echo "   featurebranches= \"${featurebranches}\""
echo "   logdir         = ${logdir}"
echo "   logfile        = ${logfile}"
echo "   machine_name_long = $machine_name_long"
echo "   prdash         = $prdash"
echo "   projects       = \"${projects}\""
echo "   regdir         = ${regdir}"
echo "   regress_mode   = ${regress_mode}"
echo "   rscriptdir     = ${rscriptdir}"
echo "   scratchdir     = ${scratchdir}"
echo " "
echo "Descriptions:"
echo "   rscriptdir -  the location of the draco regression scripts."
echo "   logdir     -  the location of the output logs."
echo "   regdir     -  the location of the top level regression system."
echo " "

# use forking to reduce total wallclock runtime, but do not fork when there is a
# dependency:
#
# draco --> capsaicin  --\
#       --> jayenne     --+--> asterisk

##---------------------------------------------------------------------------##
## Launch the jobs...
##---------------------------------------------------------------------------##

# convert featurebranches into an array
export fb=(${featurebranches})
ifb=0

# The job launch logic spawns a job for each project immediately, but the
# *-job-launch.sh script will spin until all dependencies (jobids) are met.
# Thus, the cts1-job-launch.sh for milagro will start immediately, but it will not
# do any real work until both draco and clubimc have completed.

# More sanity checks
if ! [[ -x ${rscriptdir}/${machine_class}-job-launch.sh ]]; then
   echo "FATAL ERROR: I cannot find ${rscriptdir}/${machine_class}-job-launch.sh."
   exit 1
fi

export subproj=draco
if [[ `echo $projects | grep -c $subproj` -gt 0 ]]; then
  export featurebranch=${fb[$ifb]}
  cmd="${rscriptdir}/${machine_class}-job-launch.sh"
  cmd+=" &> ${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params_sort_safe}${prdash}${featurebranch}-joblaunch.log"
  echo "${subproj}: $cmd"
  eval "${cmd} &"
  sleep 1s
  draco_jobid=`jobs -p | sort -gr | head -n 1`
  ((ifb++))
fi

export subproj=jayenne
if [[ `echo $projects | grep -c $subproj` -gt 0 ]]; then
  export featurebranch=${fb[$ifb]}
  # Run the *-job-launch.sh script (special for each platform).
  cmd="${rscriptdir}/${machine_class}-job-launch.sh"
  # Spin until $draco_jobid disappears (indicates that draco has been
  # built and installed)
  cmd+=" ${draco_jobid}"
  # Log all output.
  cmd+=" &> ${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params_sort_safe}${prdash}${featurebranch}-joblaunch.log"
  echo "${subproj}: $cmd"
  eval "${cmd} &"
  sleep 1s
  jayenne_jobid=`jobs -p | sort -gr | head -n 1`
  ((ifb++))
fi

export subproj=core
if [[ `echo $projects | grep -c $subproj` -gt 0 ]]; then
  export featurebranch=${fb[$ifb]}
  cmd="${rscriptdir}/${machine_class}-job-launch.sh"
  # Wait for draco regressions to finish
  case $extra_params_sort_safe in
  *coverage*)
     # We can only run one instance of bullseye at a time - so wait
     # for capsaicin to finish before starting core.
     cmd+=" ${draco_jobid} ${capsaicin_jobid}" ;;
  *)
     cmd+=" ${draco_jobid}" ;;
  esac
  cmd+=" &> ${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params_sort_safe}${prdash}${featurebranch}-joblaunch.log"
  echo "${subproj}: $cmd"
  eval "${cmd} &"
  sleep 1s
  core_jobid=`jobs -p | sort -gr | head -n 1`
  ((ifb++))
fi

export subproj=trt
if [[ `echo $projects | grep -c $subproj` -gt 0 ]]; then
  export featurebranch=${fb[$ifb]}
  cmd="${rscriptdir}/${machine_class}-job-launch.sh"
  # Wait for draco regressions to finish
  case $extra_params_sort_safe in
  *coverage*)
     # We can only run one instance of bullseye at a time - so wait
     # for capsaicin to finish before starting core.
     cmd+=" ${draco_jobid} ${core_jobid}" ;;
  *)
     cmd+=" ${draco_jobid} ${core_jobid}" ;;  ## trt --> core --> draco
  esac
  cmd+=" &> ${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params_sort_safe}${prdash}${featurebranch}-joblaunch.log"
  echo "${subproj}: $cmd"
  eval "${cmd} &"
  sleep 1s
  core_jobid=`jobs -p | sort -gr | head -n 1`
  ((ifb++))
fi

export subproj=npt
if [[ `echo $projects | grep -c $subproj` -gt 0 ]]; then
  export featurebranch=${fb[$ifb]}
  cmd="${rscriptdir}/${machine_class}-job-launch.sh"
  # Wait for draco regressions to finish
  case $extra_params_sort_safe in
  *coverage*)
     # We can only run one instance of bullseye at a time - so wait
     # for capsaicin to finish before starting core.
     cmd+=" ${draco_jobid} ${core_jobid} ${trt_jobid}" ;;
  *)
     cmd+=" ${draco_jobid} ${core_jobid}" ;;  ## npt --> core --> draco
  esac
  cmd+=" &> ${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params_sort_safe}${prdash}${featurebranch}-joblaunch.log"
  echo "${subproj}: $cmd"
  eval "${cmd} &"
  sleep 1s
  core_jobid=`jobs -p | sort -gr | head -n 1`
  ((ifb++))
fi

# Wait for all parallel jobs to finish
#while [ 1 ]; do fg 2> /dev/null; [ $? == 1 ] && break; done

# Wait for all subprocesses to finish before exiting this script
if [[ `jobs -p | wc -l` -gt 0 ]]; then
  echo " "
  echo "Jobs still running (if any):"
  for job in `jobs -p`; do
    echo "  waiting for job $job to finish..."
    wait $job
    echo "  waiting for job $job to finish...done"
  done
fi

# set permissions
chgrp -R ccsrad ${logdir} &> /dev/null
chmod -R g+rX,o-rwX ${logdir} &> /dev/null

echo " "
echo "All done"

##---------------------------------------------------------------------------##
## End of regression-master.sh
##---------------------------------------------------------------------------##
