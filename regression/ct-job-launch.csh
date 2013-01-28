#!/bin/tcsh

# called from ct-master.csh
# assumes the following variables are defined in ct-master.csh:
#    $regdir     - /usr/projects/jayenne/regress
#    $subproj    - 'draco', 'clubimc', etc.
#    $build_type - 'Debug', 'Release'

# sanity check
if( "${regdir}x" == "x" ) then
    echo "FATAL ERROR in ct-job-launch.csh: You did not set 'regdir' in the environment!"
    exit 1
endif
if( "${subproj}x" == "x" ) then
    echo "FATAL ERROR in ct-job-launch.csh: You did not set 'subproj' in the environment!"
    exit 1
endif
if( "${build_type}x" == "x" ) then
    echo "FATAL ERROR in ct-job-launch.csh: You did not set 'build_type' in the environment!"
    exit 1
endif

# Banner
echo "==========================================================================="
echo "CT regression job launcher for ${subproj} - ${build_type} flavor."
echo "==========================================================================="

# Configure, Build on front end
echo " "
echo "Configure and Build on the front end..."
echo "${regdir}/draco/regression/ct-regress.msub "
echo "    >& ${regdir}/logs/ct-${build_type}-${subproj}-${extra_params}-cb.log"
${regdir}/draco/regression/ct-regress.msub \
>& ${regdir}/logs/ct-${build_type}-${subproj}-${extra_params}-cb.log

# Wait for CB (Configure and Build) before starting the testing and
# reporting from the login node:

echo " "
echo "Test and Submit from the login node..."
echo "/opt/MOAB/default/bin/msub -j oe -V "
echo "    -o ${regdir}/logs/ct-${build_type}-${subproj}-${extra_params}-ts.log"
echo "    ${regdir}/draco/regression/ct-regress.msub"
set jobid = `/opt/MOAB/default/bin/msub -j oe -V -o ${regdir}/logs/ct-${build_type}-${subproj}-${extra_params}-ts.log ${regdir}/draco/regression/ct-regress.msub`

# sleep 1m
# while( `showq | grep $jobid` != "" )
#    showq | grep $jobid
#    sleep 10m
# end

# Submit from the front end
echo "Jobs done on ct."



