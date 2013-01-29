#!/bin/bash

# Use:
# - Call from crontab using
#   <path>/ccscs-debug-master.sh <build_type>
#
# - <build_type>   = {Debug, Release}
# - <extra_params> = coverage

##---------------------------------------------------------------------------##
## Support functions
##---------------------------------------------------------------------------##
print_error()
{
    echo "FATAL ERROR"
    echo " "
    echo "Usage: $0 <build_type> [extra_params]"
    echo " "
    echo "   <build_type>   = { Debug, Release }."
    echo "   [extra_params] = { coverage }."
    exit 1
}

##---------------------------------------------------------------------------##
## Main
##---------------------------------------------------------------------------##

#Banner
echo "==========================================================================="
echo "ccscs-master.sh: Regression for Linux64 on CCS LAN"

# Arguments
case $# in
1 )
    export build_type=$1
    export extra_params=none
    echo "Build: ${build_type}"
;;
2 )
    export build_type=$1
    export extra_params=$2
    echo "Build: ${build_type}     Extra Params: $extra_params"
;;
* )
    echo "Wrong number of arguments provided to ccscs-master.csh."
    print_error
    ;; 
esac

date
echo "==========================================================================="
echo " "

export regdir=/home/regress
export host=`uname -n | sed -e 's/[.].*//g'`

case ${host} in
ccscs[0-9])
    # no-op
    ;;
*)
    echo "I don't know how to run regression on host = ${host}."
    print_error
    ;;
esac

projects=(  "draco" "capsaicin" "clubimc" "wedgehog" "milagro" )
forkbuild=( "no"    "yes"       "no"      "yes"      "yes" )


# use forking to reduce total wallclock runtime, but do not fork
# when there is a dependency:
# 
# draco --> capsaicin
#       --> clubimc --> wedgehog
#                   --> milagro

# special cases
case $extra_params in
coverage)
    # no-op right now.
    ;;
# intel13)
#     # RNG fails for Release+Intel-13 so no Jayenne software
#     if ( "${build_type}" == "Release" ) then
#         set projects = (  "draco" "capsaicin" )
#         set forkbuild = ( "no"    "yes"       )
#     endif
#     breaksw
# case pgi)
#     # Capsaicin does not support building with PGI (lacking vendor installations!)
#     set projects = (  "draco" "clubimc" "wedgehog" "milagro" )
#     set forkbuild = ( "no"    "no"      "yes"      "yes" )
#     breaksw
*)
    # No-op
    ;;
esac

for (( i=0 ; i < ${#projects[@]} ; ++i )); do

    export subproj=${projects[$i]}
    export fork=${forkbuild[$i]}

    # ccscs-job-launch.csh requires the following variables:
    # $regdir, $subproj, $build_type

    echo " "
    echo "Regression for ${subproj} (${build_type}, fork=${fork})."
    echo " "
    echo "${regdir}/draco/regression/ccscs-job-launch.sh "
    echo ">& ${regdir}/logs/ccscs-${build_type}-${subproj}-${extra_params}-joblaunch.log"
        
    if test $fork = "yes"; then
        ${regdir}/draco/regression/ccscs-job-launch.sh \
        >& ${regdir}/logs/ccscs-${build_type}-${subproj}-${extra_params}-joblaunch.log &
    else
        ${regdir}/draco/regression/ccscs-job-launch.sh \
        >& ${regdir}/logs/ccscs-${build_type}-${subproj}-${extra_params}-joblaunch.log
    fi
done

##---------------------------------------------------------------------------##
## End of ccscs-master.sh
##---------------------------------------------------------------------------##
