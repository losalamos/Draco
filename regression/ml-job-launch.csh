#!/bin/tcsh

# called from ml-master.csh
# assumes the following variables are defined in ml-master.csh:
#    $regdir     - /usr/projects/jayenne/regress
#    $subproj    - 'draco', 'clubimc', etc.
#    $build_type - 'Debug', 'Release'

# sanity check
if( "${regdir}x" == "x" ) then
    echo "FATAL ERROR in ml-job-launch.csh: You did not set 'regdir' in the environment!"
    exit 1
endif
if( "${subproj}x" == "x" ) then
    echo "FATAL ERROR in ml-job-launch.csh: You did not set 'subproj' in the environment!"
    exit 1
endif
if( "${build_type}x" == "x" ) then
    echo "FATAL ERROR in ml-job-launch.csh: You did not set 'build_type' in the environment!"
    exit 1
endif

# Banner
echo "==========================================================================="
echo "ML regression job launcher for ${subproj} - ${build_type} flavor."
echo "==========================================================================="

# Configure, Build, Test on back end
echo "/opt/MOAB/bin/msub -A access -j oe -V "
echo "    -o ${regdir}/logs/ml-${build_type}-${subproj}-${extra_params}-cbt.log "
echo "    ${regdir}/draco/regression/ml-regress.msub"
set jobid = `/opt/MOAB/bin/msub -A access -j oe -V -o ${regdir}/logs/ml-${build_type}-${subproj}-${extra_params}-cbt.log ${regdir}/draco/regression/ml-regress.msub`

# Wait for CBT (Config, build, test) to finish
sleep 1m
while( `showq | grep $jobid` != "" )
   showq | grep $jobid
   sleep 10m
end

# Submit from the front end
echo "Jobs done, now submitting ${build_type} results from ml."
echo "${regdir}/draco/regression/ml-regress.msub "
echo "    >& ${regdir}/logs/ml-${build_type}-${subproj}-${extra_params}-s.log"
${regdir}/draco/regression/ml-regress.msub \
>& ${regdir}/logs/ml-${build_type}-${subproj}-${extra_params}-s.log


