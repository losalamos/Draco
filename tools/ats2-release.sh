#!/bin/bash -l

# export dry_run=1
echo -e "\nConfigure, build, and test $package for $buildflavor-$rttversion."
echo
export steps="config build test"
export logfile="$source_prefix/logs/release-$buildflavor-$rttversion-cbt${jobnameext}.log"
cmd="bsub $job_launch_options -nnodes 1 -W 4:00 \
-J rel_${package:0:3}-${buildflavor: -12}-$rttversion -o $logfile -e $logfile"
echo "$cmd $draco_script_dir/release.msub"
jobid=`$cmd $draco_script_dir/release.msub`
sleep 1m
# trim extra whitespace from number
jobid=`echo ${jobid//[^0-9]/}`
export jobids="$jobid $jobids"

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
