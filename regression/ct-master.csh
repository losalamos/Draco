#!/bin/tcsh

# Use:
# - Call from crontab using
#   <path>/ct-master.csh <build_type>
#
# - <build_type> = {Debug, Release}

#Banner
echo "==========================================================================="
echo "ct-master.csh: Regression for Cielito"

# Arguments
switch ( $#argv )
case 1:
    setenv build_type $1
    setenv extra_params none
    echo "Build: ${build_type}"
    breaksw
case 2:
    setenv build_type $1
    setenv extra_params $2
    echo "Build: ${build_type}     Extra Params: $extra_params"
    breaksw
default:
    echo "Wrong number of arguments provided to ct-master.csh."
    goto error
endsw

date
echo "==========================================================================="

switch ( "${build_type}" )
case Debug:
case Release:
    # known $build_type, continue
    breaksw
default:
    echo "unsupported build_type = ${build_type}"
    goto error
endsw

# printenv

# setenv dashboard_type Experimental
setenv regdir /usr/projects/jayenne/regress
setenv host   `echo $HOST | sed -e 's/[.].*//g'`

echo " "
switch ( "${host}" )
case ct-*:
    # module purge
    # module list
    breaksw
default:
    echo "I don't know how to run regression on host = ${host}."
    goto error
    breaksw
endsw

# use forking to reduce total wallclock runtime, but do not fork
# when there is a dependency:
# 
# draco --> capsaicin
#       --> clubimc --> wedgehog
#                   --> milagro

set projects = (  "draco" "capsaicin" "clubimc" ) # "wedgehog" "milagro" )
set forkbuild = ( "no"    "yes"       "no"      ) # "yes"      "yes" )

# special cases
switch( $extra_params )
case intel13:
    # RNG fails for Release+Intel-13 so no Jayenne software
    # if ( "${build_type}" == "Release" ) then
    #     set projects = (  "draco" "capsaicin" )
    #     set forkbuild = ( "no"    "yes"       )
    # endif
    breaksw
case pgi:
    # Capsaicin does not support building with PGI (lacking vendor installations!)
    set projects = (  "draco" "clubimc" "wedgehog" "milagro" )
    set forkbuild = ( "no"    "no"      "yes"      "yes" )
    breaksw
case none;
    # No-op
    breaksw
endsw

set i = 0
while ($i < $#projects)
    @ i = $i + 1
    setenv subproj ${projects[$i]}
    set fork    = ${forkbuild[$i]}

    # ct-job-launch.csh requires the following variables:
    # $regdir, $subproj, $build_type

    echo " "
    echo "Regression for ${subproj} (${build_type}, fork=${fork})."
    echo " "
    echo "${regdir}/draco/regression/ct-job-launch.csh "
    echo ">& ${regdir}/logs/ct-${build_type}-${subproj}-${extra_params}-joblaunch.log"
        
    if( $fork == "yes" ) then
        ${regdir}/draco/regression/ct-job-launch.csh \
        >& ${regdir}/logs/ct-${build_type}-${subproj}-${extra_params}-joblaunch.log &
    else
        ${regdir}/draco/regression/ct-job-launch.csh \
        >& ${regdir}/logs/ct-${build_type}-${subproj}-${extra_params}-joblaunch.log
    endif
end

## Labels to jump to exit OK (done) or not OK (error)
done:
    exit 0

error:  
    echo "FATAL ERROR"
    echo " "
    echo "Usage: $0 <build_type> [extra_params]"
    echo " "
    echo "   <build_type>   = { Debug, Release }."
    echo "   [extra_params] = { intel13, pgi }."
    exit 1

