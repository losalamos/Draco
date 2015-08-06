#!/bin/bash

# called from regression-master.sh
# assumes the following variables are defined in regression-master.sh
#    $regdir     - /home/regress
#    $subproj    - 'draco', 'jayenne', etc.
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

export SHOWQ=/opt/MOAB/bin/showq

 # Dependencies: wait for these jobs to finish
dep_jobids=""
for (( i=0; i < $nargs ; ++i )); do
   dep_jobids="${dep_jobids} ${args[$i]} "
done

# sanity check
if test "${regdir}x" = "x"; then
    echo "FATAL ERROR in tt-job-launch.sh: You did not set 'regdir' in the environment!"
    exit 1
fi
if test "${subproj}x" = "x"; then
    echo "FATAL ERROR in tt-job-launch.sh: You did not set 'subproj' in the environment!"
    exit 1
fi
if test "${build_type}x" = "x"; then
    echo "FATAL ERROR in tt-job-launch.sh: You did not set 'build_type' in the environment!"
    exit 1
fi
if test "${logdir}x" = "x"; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'logdir' in the environment!"
    exit 1
fi

# Banner
echo "==========================================================================="
echo "Trinitite Regression job launcher for ${subproj} - ${build_type} flavor."
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
       echo "   ${subproj}: waiting for jobid = $jobid to finish."
       sleep 5m
    done
done

# Configure, Build on front end
export mode=cb
echo " "
echo "Configure and Build on the front end..."
cmd="${regdir}/draco/regression/tt-regress.msub >& ${logdir}/tt-${build_type}-${extra_params}${epdash}${subproj}-${mode}.log"
echo "${cmd}"
eval "${cmd}"

# Wait for CB (Configure and Build) before starting the testing and
# reporting from the login node:

echo " "
echo "Test from the login node..."
cmd="/opt/MOAB/bin/msub -j oe -V -o ${logdir}/tt-${build_type}-${extra_params}${epdash}${subproj}-t.log ${regdir}/draco/regression/tt-regress.msub"
echo "${cmd}"
jobid=`eval ${cmd}`
jobid=`echo $jobid | sed '/^$/d'`
echo "jobid = ${jobid}"

# Wait for testing to finish
sleep 1m
while test "`${SHOWQ} | grep $jobid`" != ""; do
   ${SHOWQ} | grep $jobid
   sleep 1m
done

# Submit from the front end
mode=s
echo "Jobs done, now submitting ${build_type} results from tt-fey."
cmd="${regdir}/draco/regression/tt-regress.msub >& ${logdir}/tt-${build_type}-${extra_params}${epdash}${subproj}-${mode}.log"
echo "${cmd}"
eval "${cmd}"

# Submit from the front end
echo "Jobs done."

##---------------------------------------------------------------------------##
## End of tt-job-launch.sh
##---------------------------------------------------------------------------##
