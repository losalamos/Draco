#!/bin/tcsh

# called from ccscs-master.csh
# assumes the following variables are defined in ccscs-master.csh:
#    $regdir     - /home/regress
#    $subproj    - 'draco', 'clubimc', etc.
#    $build_type - 'Debug', 'Release', 'Coverage'

# sanity check
if( "${regdir}x" == "x" ) then
    echo "FATAL ERROR in ccscs-job-launch.csh: You did not set 'regdir' in the environment!"
    exit 1
endif
if( "${subproj}x" == "x" ) then
    echo "FATAL ERROR in ccscs-job-launch.csh: You did not set 'subproj' in the environment!"
    exit 1
endif
if( "${build_type}x" == "x" ) then
    echo "FATAL ERROR in ccscs-job-launch.csh: You did not set 'build_type' in the environment!"
    exit 1
endif

# Banner
echo "==========================================================================="
echo "CCSCS Regression job launcher for ${subproj} - ${build_type} flavor."
echo "==========================================================================="

# Configure, Build, Test and Submit (no Torque batch system here).
echo "${regdir}/draco/regression/ccscs-regress.msub"
echo "    >& ${regdir}/logs/ccscs-${build_type}-${subproj}-${extra_params}-cbts.log"
${regdir}/draco/regression/ccscs-regress.msub \
>& ${regdir}/logs/ccscs-${build_type}-${subproj}-${extra_params}-cbts.log


