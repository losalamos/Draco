#!/bin/bash

# called from regression-master.sh
# assumes the following variables are defined in ccscs-master.sh:
#    $regdir     - /home/regress
#    $subproj    - 'draco', 'clubimc', etc.
#    $build_type - 'Debug', 'Release'
#    $extra_params - '', 'intel13', 'pgi', 'coverage'

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

# Banner
echo "==========================================================================="
echo "CCSCS Regression job launcher for ${subproj} - ${build_type} flavor."
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

# Configure, Build, Test and Submit (no Torque batch system here).
cmd="${regdir}/draco/regression/ccscs-regress.msub >& ${regdir}/logs/ccscs-${build_type}-${extra_params}${epdash}${subproj}-cbts.log"
echo "${cmd}"
eval "${cmd}"

##---------------------------------------------------------------------------##
## End of ccscs-job-launch.sh
##---------------------------------------------------------------------------##

