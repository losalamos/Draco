#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/tt-job-launch.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# called from regression-master.sh
# assumes the following variables are defined in regression-master.sh:
#    $regdir     - /scratch/regress
#    $rscriptdir - /scratch/regress/draco/regression (actually, the location
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

# import some bash functions
export rscriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $rscriptdir/scripts/common.sh

export SHOWQ=`which squeue`
export MSUB=`which sbatch`

 # Dependencies: wait for these jobs to finish
dep_jobids=""
for (( i=0; i < $nargs ; ++i )); do
   dep_jobids="${dep_jobids} ${args[$i]} "
done

# sanity checks
job_launch_sanity_checks

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

# Select haswell or knl partition
# Optional: Use -C quad,flat to select KNL mode
# sinfo -o "%45n %30b %65f" | cut -b 47-120 | sort | uniq -c
case $extra_params in
# '--gres=' - See KT's email 2017-05-30 from Paul Peltz; Hal Marshall wants a
#             different solution (see email from same day).
#knl) partition_options="-p knl -N 1 -t 8:00:00 --gres=craynetwork:0" ;;
#*)   partition_options="-N 1 -t 8:00:00 --gres=craynetwork:0" ;;
knl) partition_options="-p knl -N 1 -t 8:00:00" ;;
*)   partition_options="-N 1 -t 8:00:00" ;;
esac

# Configure, Build on front end
echo "Configure and Build on the front end..."
export REGRESSION_PHASE=cb
echo " "
cmd="${rscriptdir}/tt-regress.msub >& ${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params}${prdash}${featurebranch}-${REGRESSION_PHASE}.log"
echo "${cmd}"
eval "${cmd}"

# Wait for CB (Configure and Build) before starting the testing and
# reporting from the worker node:
echo " "
export REGRESSION_PHASE=t
echo "Test from the login node..."
echo " "
logfile=${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params}${prdash}${featurebranch}-${REGRESSION_PHASE}.log
cmd="$MSUB -o ${logfile} -J ${subproj:0:5}-${featurebranch} --exclusive ${partition_options} ${rscriptdir}/tt-regress.msub"
echo "${cmd}"
jobid=`eval ${cmd}`
# delete blank lines
#jobid=`echo $jobid | sed '/^$/d'`
# only keep the job number
jobid=`echo $jobid | sed -e 's/.*[ ]//'`
echo "jobid = ${jobid}"

# Wait for testing to finish
sleep 1m
while test "`${SHOWQ} | grep $jobid`" != ""; do
   ${SHOWQ} | grep $jobid
   sleep 5m
done

# Submit from the front end
echo " "
echo "Submit:"
export REGRESSION_PHASE=s
echo "- jobs done, now submitting ${build_type} results from tt-fey."
cmd="${rscriptdir}/tt-regress.msub >& ${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params}${prdash}${featurebranch}-${REGRESSION_PHASE}.log"
echo "${cmd}"
eval "${cmd}"

# Submit from the front end
echo "Jobs done."

##---------------------------------------------------------------------------##
## End of tt-job-launch.sh
##---------------------------------------------------------------------------##
