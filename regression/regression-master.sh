#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/regression-master.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# Use:
# - Call from crontab using
#   <path>/regression-master.sh [options]

##---------------------------------------------------------------------------##
## Environment
##---------------------------------------------------------------------------##

# Enable job control
set -m

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
export host=`uname -n | sed -e 's/[.].*//g'`

##---------------------------------------------------------------------------##
## Support functions
##---------------------------------------------------------------------------##
ccs_extra_params="belosmods bounds_checking clang coverage fulldiagnostics gcc530 gcc610 nr perfbench valgrind "
darwin_extra_params="cuda fulldiagnostics nr perfbench valgrind"
ml_extra_params="fulldiagnostics nr perfbench pgi valgrind"
sn_extra_params="fulldiagnostics gcc610 newtools nr perfbench"
tt_extra_params="fulldiagnostics knl nr perfbench"
all_extra_params=`echo $ml_extra_params $tt_extra_params $sn_extra_params $ccs_extra_params $darwin_extra_params | xargs -n1 | sort -u | xargs`

print_use()
{
  case $host in
    ccscs*) platform_extra_params=ccs_extra_params ;;
    darwin*) platform_extra_params=darwin_extra_params ;;
    ml-*) platform_extra_params=ml_extra_params ;;
    sn-*) platform_extra_params=sn_extra_params ;;
    tt-*) platform_extra_params=tt_extra_params ;;
  esac
  platform_extra_params=`echo $platform_extra_params | sed -e 's/ / | /g'`

  echo " "
  echo "Usage: ${0##*/} -b [Release|Debug] -d [Experimental|Nightly|Continuous]"
  echo "       -h -p [\"draco jayenne capsaicin\"] -r"
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
  echo "   -p    project names  = { draco, jayenne, capsaicin }"
  echo "                          This is a space delimited list within double quotes."
  echo "   -e    extra params   = { none | $platform_extra_params }"
  echo " "
  echo "Example:"
  echo "./regression-master.sh -b Release -d Nightly -p \"draco jayenne capsaicin\""
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
   draco | jayenne | capsaicin ) # known projects, continue
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

##---------------------------------------------------------------------------##
## Main
##---------------------------------------------------------------------------##

case ${host} in
ml-*)
    export machine_name_long=Moonlight
    export machine_name_short=ml
    if [[ `fn_exists module` == 0 ]]; then
      source /usr/share/Modules/init/bash
    fi
    module purge
    # Argument checks
    if [[ ${extra_params} ]]; then
        case $extra_params in
        none)  extra_params=""; epdash="" ;;
        fulldiagnostics | nr | perfbench | pgi | valgrind ) # known, continue
        ;;
        *) echo "" ;echo "FATAL ERROR: unknown extra params (-e) = ${extra_params}"
           print_use; exit 1 ;;
        esac
    fi
    ;;
sn-*)
    export machine_name_long=Snow
    export machine_name_short=sn
    if [[ `fn_exists module` == 0 ]]; then
      source /usr/share/lmod/lmod/init/bash
    fi
    module purge
    # Argument checks
    if [[ ${extra_params} ]]; then
        case $extra_params in
        none)  extra_params=""; epdash="" ;;
        fulldiagnostics | gcc610 | newtools | nr | perfbench ) # known, continue
        ;;
        *) echo "" ;echo "FATAL ERROR: unknown extra params (-e) = ${extra_params}"
           print_use; exit 1 ;;
        esac
    fi
    ;;
tt-*)
    export machine_name_long=Trinitite
    export machine_name_short=tt
    # Argument checks
    if [[ ${extra_params} ]]; then
        case $extra_params in
        none) extra_params=""; epdash="" ;;
        fulldiagnostics | knl | nr | perfbench ) # known, continue
        ;;
        *)  echo "" ;echo "FATAL ERROR: unknown extra params (-e) = ${extra_params}"
            print_use; exit 1 ;;
        esac
    fi
    ;;
ccscs[0-9])
    export machine_name_long="Linux64 on CCS LAN"
    export machine_name_short=ccscs
    # Argument checks
    if [[ ${extra_params} ]]; then
        case $extra_params in
        none)  extra_params=""; epdash="" ;;
        belosmods | bounds_checking | coverage | fulldiagnostics | nr | perfbench | valgrind ) # known, continue
        ;;
        gcc530 | clang | gcc610 ) # known, continue
        ;;
        *) echo "" ;echo "FATAL ERROR: unknown extra params (-e) = ${extra_params}"
           print_use; exit 1 ;;
        esac
    fi
    ;;
darwin*)
    export machine_name_long="Linux64 on CCS Darwin cluster"
    export machine_name_short=darwin
    # Argument checks
    if [[ ${extra_params} ]]; then
        case $extra_params in
        none)  extra_params=""; epdash="" ;;
        cuda | fulldiagnostics | nr | perfbench | valgrind ) # known, continue
        ;;
        *) echo "" ;echo "FATAL ERROR: unknown extra params (-e) = ${extra_params}"
           print_use; exit 1 ;;
        esac
    fi
    ;;
*)
    echo "FATAL ERROR: I don't know how to run regression on host = ${host}."
    print_use;  exit 1 ;;
esac

#------------------------------------------------------------------------------#
# Redirect output to a logfile.
#------------------------------------------------------------------------------#

# Ensure log dir exists.
mkdir -p $logdir || die "Could not create a directory for log files."

# Redirect output to logfile.
timestamp=`date +%Y%m%d-%H%M`
logfile=$logdir/${machine_name_short}-${build_type}-master-$timestamp.log
#echo "Redirecting output to $logfile"
exec > $logfile
exec 2>&1

##---------------------------------------------------------------------------##
## Export environment
##---------------------------------------------------------------------------##
export build_autodoc build_type dashboard_type epdash extra_params
export featurebranches logdir prdash regdir regress_mode scratchdir

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
echo "   extra_params   = ${extra_params}"
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
# Thus, the ml-job-launch.sh for milagro will start immediately, but it will not
# do any real work until both draco and clubimc have completed.

# More sanity checks
if ! [[ -x ${rscriptdir}/${machine_name_short}-job-launch.sh ]]; then
   echo "FATAL ERROR: I cannot find ${rscriptdir}/${machine_name_short}-job-launch.sh."
   exit 1
fi

export subproj=draco
if [[ `echo $projects | grep -c $subproj` -gt 0 ]]; then
  export featurebranch=${fb[$ifb]}
  cmd="${rscriptdir}/${machine_name_short}-job-launch.sh"
  cmd+=" &> ${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params}${prdash}${featurebranch}-joblaunch.log"
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
  cmd="${rscriptdir}/${machine_name_short}-job-launch.sh"
  # Spin until $draco_jobid disappears (indicates that draco has been
  # built and installed)
  cmd+=" ${draco_jobid}"
  # Log all output.
  cmd+=" &> ${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params}${prdash}${featurebranch}-joblaunch.log"
  echo "${subproj}: $cmd"
  eval "${cmd} &"
  sleep 1s
  jayenne_jobid=`jobs -p | sort -gr | head -n 1`
  ((ifb++))
fi

export subproj=capsaicin
if [[ `echo $projects | grep -c $subproj` -gt 0 ]]; then
  export featurebranch=${fb[$ifb]}
  cmd="${rscriptdir}/${machine_name_short}-job-launch.sh"
  # Wait for draco regressions to finish
  case $extra_params in
  coverage)
     # We can only run one instance of bullseye at a time - so wait
     # for jayenne to finish before starting capsaicin.
     cmd+=" ${draco_jobid} ${jayenne_jobid}" ;;
  *)
     cmd+=" ${draco_jobid}" ;;
  esac
  cmd+=" &> ${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params}${prdash}${featurebranch}-joblaunch.log"
  echo "${subproj}: $cmd"
  eval "${cmd} &"
  sleep 1s
  capsaicin_jobid=`jobs -p | sort -gr | head -n 1`
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
chgrp -R draco ${logdir} &> /dev/null
chmod -R g+rX ${logdir} &> /dev/null

echo " "
echo "All done"

##---------------------------------------------------------------------------##
## End of regression-master.sh
##---------------------------------------------------------------------------##
