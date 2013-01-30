#!/bin/bash

# Use:
# - Call from crontab using
#   <path>/regression-master.sh <build_type> <extra_params>

##---------------------------------------------------------------------------##
## Support functions
##---------------------------------------------------------------------------##
print_use()
{
    echo " "
    echo "Usage: $0 <build_type> [extra_params]"
    echo " "
    echo "   <build_type>   = { Debug, Release }."
    echo "   [extra_params] = { intel13, pgi, coverage }."
    echo " "
    echo "Extra parameters read from environment:"
    echo "   ENV{dashboard_type} = {Nightly, Experimental}"
    echo "   ENV{base_dir}       = {/var/tmp/$USER/cdash, /scratch/$USER/cdash}"
}

##---------------------------------------------------------------------------##
## Main
##---------------------------------------------------------------------------##

# Defaults
if test "${dashboard_type}x" = "x"; then
    export dashboard_type=Nightly
fi

# Arguments
case $# in
1 )
    export build_type=$1
    export extra_params=""
;;
2 )
    export build_type=$1
    export extra_params=$2
;;
* )
    echo "FATAL ERROR: Wrong number of arguments provided to regression-master.sh."
    print_use
    exit 1
    ;; 
esac

# Host based variables
export host=`uname -n | sed -e 's/[.].*//g'`

case ${host} in
ct-*)
    machine_name_long=Cielito
    machine_name_short=ct
    export regdir=/usr/projects/jayenne/regress
    # We don't include capsaicin in the Intel-12 based builds.
    projects=(  "draco" "clubimc" "wedgehog" "milagro" )
    forkbuild=( "no"    "no"      "yes"      "yes" )
    ;;
ml-*)
    machine_name_long=Moonlight
    machine_name_short=ml
    module purge
    export regdir=/usr/projects/jayenne/regress
    # We don't include capsaicin in the Intel-12 based builds.
    projects=(  "draco" "clubimc" "wedgehog" "milagro" )
    forkbuild=( "no"    "no"      "yes"      "yes" )
    ;;
ccscs[0-9])
    machine_name_long="Linux64 on CCS LAN"
    machine_name_short=ccscs
    export regdir=/home/regress
    projects=(  "draco" "capsaicin" "clubimc" "wedgehog" "milagro" )
    forkbuild=( "no"    "yes"       "no"      "yes"      "yes" )
    ;;
*)
    echo "FATAL ERROR: I don't know how to run regression on host = ${host}."
    print_use
    exit 1
    ;;
esac

# Banner

echo "==========================================================================="
echo "regression-master.sh: Regression for $machine_name_long"
echo "Build: ${build_type}     Extra Params: $extra_params"
date
echo "==========================================================================="
echo " "
echo "Environment:"
echo "   build_type   = ${build_type}"
echo "   extra_params = ${extra_params}"
echo "   regdir       = ${regdir}"
echo " "
echo "Optional environment:"
echo "   dashboard_type = ${dashboard_type}"
echo "   base_dir       = ${base_dir}"
echo " "

# Sanity Checks
case ${build_type} in
"Debug" | "Release" )
    # known $build_type, continue
    ;;
*)
    echo "FATAL ERROR: unsupported build_type = ${build_type}"
    print_use
    exit 1
    ;; 
esac

case ${dashboard_type} in
Nightly | Experimental)
    # known dashboard_type, continue
    ;;
*)
    echo "FATAL ERROR: unknown dashboard_type = ${dashboard_type}"
    print_use
    exit 1
    ;;
esac

# use forking to reduce total wallclock runtime, but do not fork
# when there is a dependency:
# 
# draco --> capsaicin
#       --> clubimc --> wedgehog
#                   --> milagro

# special cases
case $extra_params in
intel13)
    # also build capsaicin
    projects=(  "draco" "capsaicin" "clubimc" "wedgehog" "milagro" )
    forkbuild=( "no"    "yes"       "no"      "yes"      "yes" )
    epdash="-"
    ;;
pgi)
    # Capsaicin does not support building with PGI (lacking vendor installations!)
    projects=(  "draco" "clubimc" "wedgehog" "milagro" )
    forkbuild=( "no"    "no"      "yes"      "yes" )
    epdash="-"
    ;;
coverage)
    epdash="-"
    ;;
*)
    epdash=""
    ;;
esac

for (( i=0 ; i < ${#projects[@]} ; ++i )); do

    export subproj=${projects[$i]}
    export fork=${forkbuild[$i]}

    # <machine>-job-launch.sh requires the following variables:
    # $regdir, $subproj, $build_type

    cmd="${regdir}/draco/regression/${machine_name_short}-job-launch.sh >& ${regdir}/logs/${machine_name_short}-${build_type}-${extra_params}${epdash}${subproj}-joblaunch.log"

    echo " "
    echo "Regression for ${subproj} (${build_type}, fork=${fork})."
    echo " "
    echo "${cmd}"
        
    if test $fork = "yes"; then
        eval ${cmd} &
    else
        eval ${cmd}
    fi
done

##---------------------------------------------------------------------------##
## End of regression-master.sh
##---------------------------------------------------------------------------##
