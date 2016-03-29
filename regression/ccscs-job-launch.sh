#!/bin/bash

# called from regression-master.sh
# assumes the following variables are defined in regression-master.sh:
#    $regdir     - /home/regress
#    $subproj    - 'draco', 'clubimc', 'jaynne', etc.
#    $build_type - 'Debug', 'Release'
#    $extra_params - '', 'intel13', 'pgi', 'coverage'

# command line arguments
args=( "$@" )
nargs=${#args[@]}

 # Dependencies: wait for these jobs to finish
dep_jobids=""
for (( i=0; i < $nargs ; ++i )); do
   dep_jobids="${dep_jobids} ${args[$i]} "
done

# sanity check
if test "${regdir}x" = "x"; then
    echo "FATAL ERROR in ccscs-job-launch.sh: You did not set 'regdir' in the environment!"
    exit 1
fi
if test "${subproj}x" = "x"; then
    echo "FATAL ERROR in ccscs-job-launch.sh: You did not set 'subproj' in the environment!"
    exit 1
fi
if test "${build_type}x" = "x"; then
    echo "FATAL ERROR in ccscs-job-launch.sh: You did not set 'build_type' in the environment!"
    exit 1
fi
if test "${logdir}x" = "x"; then
    echo "FATAL ERROR in ccscs-job-launch.sh: You did not set 'logdir' in the environment!"
    exit 1
fi

# Banner
echo "==========================================================================="
echo "CCSCS Regression job launcher"
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

# epdash="-"
# if test "${extra_params}x" = "x"; then
#    epdash=""
# fi

# Prerequisits:
# Wait for all dependencies to be met before creating a new job

for jobid in ${dep_jobids}; do
    while [ `ps --no-headers -u ${USER} -o pid | grep ${jobid} | wc -l` -gt 0 ]; do
       echo "   ${subproj}: waiting for jobid = $jobid to finish (sleeping 5 min)."
       sleep 5m
    done
done

# Configure, Build, Test and Submit (no Torque batch system here).
cmd="${regdir}/draco/regression/ccscs-regress.msub >& ${logdir}/ccscs-${build_type}-${extra_params}${epdash}${subproj}-cbts.log"
echo "${cmd}"
eval "${cmd}"

echo "All done."

##---------------------------------------------------------------------------##
## End of ccscs-job-launch.sh
##---------------------------------------------------------------------------##
