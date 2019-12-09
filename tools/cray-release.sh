#!/bin/bash

# export dry_run=1
# config and build on front-end
echo -e "\nConfigure and build $package for $buildflavor-$rttversion."
echo
export steps="config build"
logfile="$source_prefix/logs/release-$buildflavor-$rttversion-cb.log"
run "$draco_script_dir/release.msub &> $logfile"

# Run the tests on the back-end.
export steps="test"
logfile="$source_prefix/logs/release-$buildflavor-$rttversion-t.log"
cmd="sbatch -J rel-${package:0:3}-${buildflavor: -10}-$rttversion -t 8:00:00 -N 1 \
 --gres=craynetwork:0 -o $logfile $draco_script_dir/release.msub"
echo -e "\nTest $package for $buildflavor-$rttversion."
echo "$cmd"
jobid=`eval ${cmd}`
sleep 1m
# trim extra whitespace from number
jobid=`echo ${jobid//[^0-9]/}`
export jobids="$jobid $jobids"

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
