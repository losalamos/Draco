#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/cts1-job-launch.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2019, Triad National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# called from regression-master.sh
# assumes the following variables are defined in regression-master.sh:
#    $regdir     - /usr/projects/jayenne/regress
#    $rscriptdir - $regdir/draco/regression (actually, the location
#                  where the active regression_master.sh is located)
#    $subproj    - 'draco', 'jaynne', 'capsaicin', etc.
#    $build_type - 'Debug', 'Release'
#    $extra_params - '', 'intel13', 'pgi', 'coverage'

# Under cron, a basic environment might not be loaded yet.
if [[ `which sbatch 2>/dev/null | grep -c sbatch` == 0 ]]; then
  source /etc/bash.bashrc.local
fi

# command line arguments
args=( "$@" )
nargs=${#args[@]}
scriptname=${0##*/}
host=`uname -n`

export SHOWQ=`which squeue`
export MSUB=`which sbatch`

# Dependencies: wait for these jobs to finish
dep_jobids=""
for (( i=0; i < $nargs ; ++i )); do
   dep_jobids="${dep_jobids} ${args[$i]} "
done

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

# sanity checks
job_launch_sanity_checks

available_queues=`sacctmgr -np list assoc user=$LOGNAME | sed -e 's/.*|\(.*dev.*\)|.*/\1/' | sed -e 's/|.*//'`
available_account=`sacctmgr -np list assoc user=$LOGNAME format=Account | sed -e 's/|//' | sed -e 's/,/ /g' | xargs -n 1 | sort -u | xargs`
available_partition=`sacctmgr -np list assoc user=$LOGNAME format=Partition | sed -e 's/|//' | sed -e 's/,/ /g' | xargs -n 1 | sort -u | xargs`
available_qos=`sacctmgr -np list assoc user=$LOGNAME format=Qos | sed -e 's/|//' | sed -e 's/,/ /g' | xargs -n 1 | sort -u | xargs`
# not sure how to detect access to
# pm="--reservation=PreventMaint"
case $available_queues in
  *access*) access_queue="-A access --qos=access" ;;
  *dev*)    access_queue="--qos=dev" ;;
esac

case $machine_name_short in
  # Must use interactive qos on grizzly because standrad has a 70 node minimum
  gr*) access_queue="--qos=interactive";;
esac

# Banner
print_job_launch_banner

# Prerequisits:
# Wait for all dependencies to be met before creating a new job

for jobid in ${dep_jobids}; do
    while [ `ps --no-headers -u ${USER} -o pid | grep -c ${jobid}` -gt 0 ]; do
       echo "   ${subproj}: waiting for jobid = $jobid to finish (sleeping 5 min)."
       sleep 5m
    done
done

if ! [[ -d $logdir ]]; then
  mkdir -p $logdir
  chgrp draco $logdir
  chmod g+rwX $logdir
  chmod g+s $logdir
fi

build_partition_options="-N 1 -t 4:00:00"
partition_options="-N 1 -t 4:00:00"

# Configure on the front end
echo "Configure:"
export REGRESSION_PHASE=c
cmd="${rscriptdir}/cts1-regress.msub >& ${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params_sort_safe}${prdash}${featurebranch}-${REGRESSION_PHASE}.log"
echo "${cmd}"
eval "${cmd}"

# Build, Test on back end
echo " "
echo "Build, Test:"
export REGRESSION_PHASE=bt
logfile=${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params_sort_safe}${prdash}${featurebranch}-${REGRESSION_PHASE}.log
if [[ -f $logfile ]]; then
  rm $logfile
fi
cmd="$MSUB ${access_queue} -o ${logfile} -J ${subproj:0:5}-${featurebranch} ${build_partition_options} ${rscriptdir}/${machine_class}-regress.msub"
echo "${cmd}"
jobid=`eval ${cmd}`
# trim extra whitespace from number
jobid=`echo ${jobid//[^0-9]/}`

# Wait for BT (build and test) to finish
sleep 1m
while test "`$SHOWQ | grep $jobid`" != ""; do
   $SHOWQ | grep $jobid
   echo "   ${subproj}: waiting for jobid = $jobid to finish (sleeping 5 minutes)."
   sleep 5m
done

# Submit from the front end
echo " "
echo "Submit:"
export REGRESSION_PHASE=s
echo "Jobs done, now submitting ${build_type} results from ${host}."
cmd="${rscriptdir}/cts1-regress.msub >& ${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params_sort_safe}${prdash}${featurebranch}-s.log"
echo "${cmd}"
eval "${cmd}"

echo "All done."

##---------------------------------------------------------------------------##
## End of ba-job-launch.sh
##---------------------------------------------------------------------------##
