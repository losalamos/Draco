#!/bin/bash

# called from ml-master.cs
# assumes the following variables are defined in ml-master.cs:
#    $regdir     - /home/regress
#    $subproj    - 'draco', 'clubimc', etc.
#    $build_type - 'Debug', 'Release', 'Coverage'

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

# Configure, Build, Test on back end
cmd="/opt/MOAB/bin/msub -A access -j oe -V -o ${regdir}/logs/ml-${build_type}-${extra_params}${epdash}${subproj}-cbt.log ${regdir}/draco/regression/ml-regress.msub"
echo "${cmd}"
jobid=`eval ${cmd}`

# Wait for CBT (Config, build, test) to finish
sleep 1m
while test "`showq | grep $jobid`" != ""; do
   showq | grep $jobid
   sleep 10m
done

# Submit from the front end
echo "Jobs done, now submitting ${build_type} results from ml."
cmd="${regdir}/draco/regression/ml-regress.msub >& ${regdir}/logs/ml-${build_type}-${extra_params}${epdash}${subproj}-s.log"
echo "${cmd}"
eval "${cmd}"

##---------------------------------------------------------------------------##
## End of ml-job-launch.sh
##---------------------------------------------------------------------------##

