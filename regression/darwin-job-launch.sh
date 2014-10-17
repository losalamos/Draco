#!/bin/bash

# called from regression-master.sh
# assumes the following variables are defined in regression-master.sh:
#    $regdir     - /home/regress
#    $subproj    - 'draco', 'clubimc', etc.
#    $build_type - 'Debug', 'Release', 'Coverage'

# command line arguments
args=( "$@" )
nargs=${#args[@]}
scriptname=`basename $0`
host=`uname -n`

# if test ${nargs} -lt 1; then
#     echo "Fatal Error: launch job requires a subproject name"
#     echo " "
#     echo "Use:"
#     echo "   launchjob projname [jobid] [jobid]"
#     echo " "
#     return 1
#     # exit 1
# fi

#export MOABHOMEDIR=/opt/MOAB
export SHOWQ=/usr/bin/squeue

# Dependencies: wait for these jobs to finish
dep_jobids=""
for (( i=0; i < $nargs ; ++i )); do
   dep_jobids="${dep_jobids} ${args[$i]} "
done

# sanity check
if test "${regdir}x" = "x"; then
    echo "FATAL ERROR in darwin-job-launch.sh: You did not set 'regdir' in the environment!"
    exit 1
fi
if test "${subproj}x" = "x"; then
    echo "FATAL ERROR in darwin-job-launch.sh: You did not set 'subproj' in the environment!"
    exit 1
fi
if test "${build_type}x" = "x"; then
    echo "FATAL ERROR in darwin-job-launch.sh: You did not set 'build_type' in the environment!"
    exit 1
fi
if test "${logdir}x" = "x"; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'logdir' in the environment!"
    exit 1
fi

# Banner
echo "==========================================================================="
echo "Darwin Regression job launcher."
echo "==========================================================================="
echo " "
echo "Environment:"
echo "   subproj        = ${subproj}"
echo "   build_type     = ${build_type}"
if test "${extra_params}x" == "x"; then
echo "   extra_params   = none"
else
echo "   extra_params   = ${extra_params}"
fi
echo "   regdir         = ${regdir}"
echo "   logdir         = ${logdir}"
echo "   dashboard_type = ${dashboard_type}"
#echo "   base_dir       = ${base_dir}"
echo " "
echo "   ${subproj}: dep_jobids = ${dep_jobids}"
echo " "

# Prerequisits:
# Wait for all dependencies to be met before creating a new job

for jobid in ${dep_jobids}; do
    while [ `ps --no-headers -u ${USER} -o pid | grep ${jobid} | wc -l` -gt 0 ]; do
       echo "   ${subproj}: waiting for jobid = $jobid to finish (sleeping 5 minutes)."
       sleep 5m
    done
done

# Configure, Build, Test on back end
cmd="/usr/bin/sbatch -v -o ${logdir}/darwin-${build_type}-${extra_params}${epdash}${subproj}-cbt.log -e ${regdir}/logs/darwin-${build_type}-${extra_params}${epdash}${subproj}-cbt.log ${regdir}/draco/regression/darwin-regress.msub"
echo "${cmd}"
jobid=`eval ${cmd}`
# trim extra whitespace from number
jobid=`echo ${jobid//[^0-9]/}`

# Wait for CBT (Config, build, test) to finish
sleep 1m
while test "`$SHOWQ | grep $jobid`" != ""; do
   $SHOWQ | grep $jobid
   echo "   ${subproj}: waiting for jobid = $jobid to finish (sleeping 5 minutes)."
   sleep 5m
done

# Submit from the front end
echo "Jobs done, now submitting ${build_type} results from darwin."
cmd="${regdir}/draco/regression/darwin-regress.msub >& ${logdir}/darwin-${build_type}-${extra_params}${epdash}${subproj}-s.log"
echo "${cmd}"
eval "${cmd}"

echo "All done."

##---------------------------------------------------------------------------##
## End of darwin-job-launch.sh
##---------------------------------------------------------------------------##

