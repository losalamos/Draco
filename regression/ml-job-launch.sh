#!/bin/bash

# called from ml-master.cs
# assumes the following variables are defined in ml-master.cs:
#    $regdir     - /home/regress
#    $subproj    - 'draco', 'clubimc', etc.
#    $build_type - 'Debug', 'Release', 'Coverage'

# command line arguments
args=( "$@" )
nargs=${#args[@]}

# if test ${nargs} -lt 1; then
#     echo "Fatal Error: launch job requires a subproject name"
#     echo " "
#     echo "Use:"
#     echo "   launchjob projname [jobid] [jobid]"
#     echo " "
#     return 1
#     # exit 1
# fi

export MOABHOMEDIR=/opt/MOAB
export SHOWQ=/opt/MOAB/bin/showq

# Dependencies: wait for these jobs to finish
dep_jobids=""
for (( i=0; i < $nargs ; ++i )); do
   dep_jobids="${dep_jobids} ${args[$i]} "
done

# sanity check
if test "${regdir}x" = "x"; then
    echo "FATAL ERROR in ml-job-launch.sh: You did not set 'regdir' in the environment!"
    exit 1
fi
if test "${subproj}x" = "x"; then
    echo "FATAL ERROR in ml-job-launch.sh: You did not set 'subproj' in the environment!"
    exit 1
fi
if test "${build_type}x" = "x"; then
    echo "FATAL ERROR in ml-job-launch.sh: You did not set 'build_type' in the environment!"
    exit 1
fi

# Banner
echo "==========================================================================="
echo "ML Regression job launcher for ${subproj} - ${build_type} flavor."
echo "==========================================================================="
echo " "
echo "Environment:"
echo "   subproj      = ${subproj}"
echo "   build_type   = ${build_type}"
echo "   extra_params = ${extra_params}"
echo "   regdir       = ${regdir}"
echo " "
echo "Optional environment:"
echo "   dashboard_type = ${dashboard_type}"
echo "   base_dir       = ${base_dir}"
echo " "

epdash="-"
if test "${extra_params}x" = "x"; then
   epdash=""
fi

# Prerequisits:
# Wait for all dependencies to be met before creating a new job
echo "   ${subproj}: dep_jobids = ${dep_jobids}"
for jobid in ${dep_jobids}; do
    while [ `ps --no-headers -u ${USER} -o pid | grep ${jobid} | wc -l` -gt 0 ]; do
       echo "   ${subproj}: waiting for jobid = $jobid to finish."
       sleep 10m
    done
done

# Configure, Build, Test on back end
cmd="/opt/MOAB/bin/msub -A access -j oe -V -o ${regdir}/logs/ml-${build_type}-${extra_params}${epdash}${subproj}-cbt.log ${regdir}/draco/regression/ml-regress.msub"
echo "${cmd}"
jobid=`eval ${cmd}`
# trim extra whitespace from number
jobid=`echo ${jobid//[^0-9]/}`

# Wait for CBT (Config, build, test) to finish
sleep 1m
while test "`$SHOWQ | grep $jobid`" != ""; do
   $SHOWQ | grep $jobid
   sleep 10m
done

# Submit from the front end
echo "Jobs done, now submitting ${build_type} results from ml."
cmd="${regdir}/draco/regression/ml-regress.msub >& ${regdir}/logs/ml-${build_type}-${extra_params}${epdash}${subproj}-s.log"
echo "${cmd}"
eval "${cmd}"

echo "All done."

##---------------------------------------------------------------------------##
## End of ml-job-launch.sh
##---------------------------------------------------------------------------##

