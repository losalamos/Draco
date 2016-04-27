#!/bin/bash

# called from regression-master.cs
# assumes the following variables are defined in regression-master.cs:
#    $regdir     - /home/regress
#    $subproj    - 'draco', 'clubimc', 'jaynne', etc.
#    $build_type - 'Debug', 'Release'
#    $extra_params - '', 'intel13', 'pgi', 'coverage'

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

export MOABHOMEDIR=/opt/MOAB
extradirs="/opt/MOAB/bin /opt/MOAB/default/bin"
for mydir in ${extradirs}; do
   if test -z "`echo $PATH | grep $mydir`" && test -d $mydir; then
      export PATH=${PATH}:${mydir}
   fi
done
export SHOWQ=`which showq`

# Dependencies: wait for these jobs to finish
dep_jobids=""
for (( i=0; i < $nargs ; ++i )); do
   dep_jobids="${dep_jobids} ${args[$i]} "
done

# sanity check
if test "${regdir}x" = "x"; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'regdir' in the environment!"
    exit 1
fi
if test "${subproj}x" = "x"; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'subproj' in the environment!"
    exit 1
fi
if test "${build_type}x" = "x"; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'build_type' in the environment!"
    exit 1
fi
if test "${logdir}x" = "x"; then
    echo "FATAL ERROR in ${scriptname}: You did not set 'logdir' in the environment!"
    exit 1
fi

# What queue should we use
#access_queue=""
#if test -x /opt/MOAB/bin/drmgroups; then
#   avail_queues=`/opt/MOAB/bin/drmgroups`
avail_queues=`mdiag -u $LOGNAME | grep ALIST | sed -e 's/.*ALIST=//' | sed -e 's/,/ /g'`
case $avail_queues in
*access*) access_queue="-A access" ;;
esac
#fi

# Banner
echo "==========================================================================="
echo "ML Regression job launcher"
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
echo "   MOAB queue     = ${access_queue}"
echo " "
echo "   ${subproj}: dep_jobids = ${dep_jobids}"
echo " "

echo "module purge &> /dev/null"
module purge &> /dev/null
echo "module list"
module list

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

# Configure, Build, Test on back end
cmd="/opt/MOAB/bin/msub ${access_queue} -j oe -V -o ${logdir}/ml-${build_type}-${extra_params}${epdash}${subproj}-cbt.log ${regdir}/draco/regression/ml-regress.msub"
echo "${cmd}"
jobid=`eval ${cmd}`
# trim extra whitespace from number
jobid=`echo ${jobid//[^0-9]/}`

# Wait for CBT (Config, build, test) to finish
sleep 1m
while test "`$SHOWQ | grep $jobid`" != ""; do
   $SHOWQ | grep $jobid
   sleep 5m
done

# Submit from the front end
echo "Jobs done, now submitting ${build_type} results from ${host}."
cmd="${regdir}/draco/regression/ml-regress.msub >& ${logdir}/ml-${build_type}-${extra_params}${epdash}${subproj}-s.log"
echo "${cmd}"
eval "${cmd}"

echo "All done."

##---------------------------------------------------------------------------##
## End of ml-job-launch.sh
##---------------------------------------------------------------------------##
