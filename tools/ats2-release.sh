#!/bin/bash -l

# export dry_run=1
echo -ne "\nConfigure, build, and test ${package:-unknown} for "
echo -ne "${buildflavor:-unknown}-${rttversion:-unknown}.\n"
export steps="config build test"
export logfile="${source_prefix:-unknoan}/logs/release-$buildflavor-$rttversion\
-cbt${jobnameext:-unknown}.log"
# shellcheck disable=SC2154
cmd="bsub ${job_launch_options} -nnodes 1 -W 8:00 \
-J rel_${package:0:3}-${buildflavor: -12}-$rttversion -o $logfile -e $logfile"
echo "$cmd ${draco_script_dir:-unknown}/release.msub"
jobid=$(${cmd} "$draco_script_dir/release.msub")
echo "jobid = $jobid"
sleep 2m
# trim extra whitespace from number
jobid="${jobid//[^0-9]/}"
export jobids="$jobid $jobids"
echo "jobids = $jobids"

# ------------------------------------------------------------------------------------------------ #
# End ats2-release.sh
# ------------------------------------------------------------------------------------------------ #
