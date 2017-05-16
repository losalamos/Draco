#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/ccscs-job-launch.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# called from regression-master.sh
# assumes the following variables are defined in regression-master.sh:
#    $regdir     - /scratch/regress
#    $rscriptdir - /scratch/regress/draco/regression (actually, the location
#                  where the active regression_master.sh is located)
#    $subproj    - 'draco', 'jaynne', 'capsaicin', etc.
#    $build_type - 'Debug', 'Release'
#    $extra_params - '', 'intel13', 'pgi', 'coverage'

# command line arguments
args=( "$@" )
nargs=${#args[@]}
scriptname=${0##*/}
host=`uname -n`

 # Dependencies: wait for these jobs to finish
dep_jobids=""
for (( i=0; i < $nargs ; ++i )); do
   dep_jobids="${dep_jobids} ${args[$i]} "
done

# sanity check
if [[ ! ${regdir} ]]; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'regdir' in the environment!"
    echo "printenv -> "
    printenv
    exit 1
fi
if [[ ! ${rscriptdir} ]]; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'rscriptdir' in the environment!"
    echo "printenv -> "
    printenv
    exit 1
fi
if [[ ! ${subproj} ]]; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'subproj' in the environment!"
    echo "printenv -> "
    printenv
    exit 1
fi
if [[ ! ${build_type} ]]; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'build_type' in the environment!"
    echo "printenv -> "
    printenv
    exit 1
fi
if [[ ! ${logdir} ]]; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'logdir' in the environment!"
    echo "printenv -> "
    printenv
    exit 1
fi

if test $subproj == draco || test $subproj == jayenne; then
  if [[ ! ${featurebranch} ]]; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'featurebranch' in the environment!"
    echo "printenv -> "
    printenv
  fi
fi

# Banner
echo "==========================================================================="
echo "CCSCS Regression job launcher for ${subproj} - ${build_type} flavor."
echo "==========================================================================="
echo " "
echo "Environment:"
if [[ ${featurebranch} ]]; then
  echo "   featurebranch  = ${featurebranch}"
fi
echo "   build_autodoc  = ${build_autodoc}"
echo "   build_type     = ${build_type}"
echo "   dashboard_type = ${dashboard_type}"
echo "   epdash         = $epdash"
if [[ ! ${extra_params} ]]; then
  echo "   extra_params   = none"
else
  echo "   extra_params   = ${extra_params}"
fi
if [[ ${featurebranch} ]]; then
  echo "   featurebranch  = ${featurebranch}"
fi
echo "   logdir         = ${logdir}"
echo "   logfile        = ${logfile}"
echo "   machine_name_long = $machine_name_long"
echo "   prdash         = $prdash"
echo "   projects       = \"${projects}\""
echo "   regdir         = ${regdir}"
echo "   regress_mode   = ${regress_mode}"
echo "   rscriptdir     = ${rscriptdir}"
echo "   scratchdir     = ${scratchdir}"
echo "   subproj        = ${subproj}"
echo " "
echo "   ${subproj}: dep_jobids = ${dep_jobids}"
echo " "

# Prerequisits:
# Wait for all dependencies to be met before creating a new job

for jobid in ${dep_jobids}; do
    while [ `ps --no-headers -u ${USER} -o pid | grep -c ${jobid}` -gt 0 ]; do
       echo "   ${subproj}: waiting for jobid = $jobid to finish (sleeping 5 min)."
       sleep 5m
    done
done

if ! test -d $logdir; then
  mkdir -p $logdir
  chgrp draco $logdir
  chmod g+rwX $logdir
  chmod g+s $logdir
fi

# Configure, Build, Test and Submit (no Torque batch system here).
# (c)onfigure, (b)uild, (t)est, (s)ubmit

echo "Configure, Build, Test:"
export REGRESSION_PHASE=cbt
cmd="${rscriptdir}/ccscs-regress.msub >& ${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params}${prdash}${featurebranch}-${REGRESSION_PHASE}.log"
echo "${cmd}"
eval "${cmd}"

echo "Submit:"
export REGRESSION_PHASE=s
cmd="${rscriptdir}/ccscs-regress.msub >& ${logdir}/${machine_name_short}-${subproj}-${build_type}${epdash}${extra_params}${prdash}${featurebranch}-${REGRESSION_PHASE}.log"
echo "${cmd}"
eval "${cmd}"

echo "All done."

##---------------------------------------------------------------------------##
## End of ccscs-job-launch.sh
##---------------------------------------------------------------------------##
