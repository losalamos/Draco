#!/bin/bash

# called from ct-master.cs
# assumes the following variables are defined in ct-master.cs:
#    $regdir     - /home/regress
#    $subproj    - 'draco', 'clubimc', etc.
#    $build_type - 'Debug', 'Release', 'Coverage'

# sanity check
if test "${regdir}x" = "x"; then
    echo "FATAL ERROR in ct-job-launch.sh: You did not set 'regdir' in the environment!"
    exit 1
fi
if test "${subproj}x" = "x"; then
    echo "FATAL ERROR in ct-job-launch.sh: You did not set 'subproj' in the environment!"
    exit 1
fi
if test "${build_type}x" = "x"; then
    echo "FATAL ERROR in ct-job-launch.sh: You did not set 'build_type' in the environment!"
    exit 1
fi

# Banner
echo "==========================================================================="
echo "CT Regression job launcher for ${subproj} - ${build_type} flavor."
echo "==========================================================================="

echo " "
echo "Environment variables used by script:"
echo "  subproj      = ${subproj}"
echo "  build_type   = ${build_type}"
echo "  extra_params = ${extra_params}"
echo "  regdir       = ${regdir}"
echo " "
echo "Optional environment:"
echo "   dashboard_type = ${dashboard_type}"
echo "   base_dir       = ${base_dir}"
echo " "

epdash="-"
if test "${extra_params}x" = "x"; then
   epdash=""
fi

# Configure, Build on front end
export mode=cb
echo " "
echo "Configure and Build on the front end..."
cmd="${regdir}/draco/regression/ct-regress.msub >& ${regdir}/logs/ct-${build_type}-${extra_params}${epdash}${subproj}-${mode}.log"
echo "${cmd}"
eval "${cmd}"

# Wait for CB (Configure and Build) before starting the testing and
# reporting from the login node:

echo " "
echo "Test from the login node..."
cmd="/opt/MOAB/default/bin/msub -j oe -V -o ${regdir}/logs/ct-${build_type}-${extra_params}${epdash}${subproj}-t.log ${regdir}/draco/regression/ct-regress.msub"
echo "${cmd}"
jobid=`eval ${cmd}`
echo "jobid = ${jobid}"

# Wait for testing to finish
sleep 1m
while test "`showq | grep $jobid`" != ""; do
   showq | grep $jobid
   sleep 10m
done

# Submit from the front end
mode=s
echo "Jobs done, now submitting ${build_type} results from ct-fe1."
cmd="${regdir}/draco/regression/ct-regress.msub >& ${regdir}/logs/ct-${build_type}-${extra_params}${epdash}${subproj}-${mode}.log"
echo "${cmd}"
eval "${cmd}"


# Submit from the front end
echo "Jobs done on ct."

##---------------------------------------------------------------------------##
## End of ct-job-launch.sh
##---------------------------------------------------------------------------##

