#!/bin/bash -l

##---------------------------------------------------------------------------##
## Assumptions:
##---------------------------------------------------------------------------##
## 1. Directory layout:
##    /usr/projects/draco/draco-NN_NN_NN/
##                  scripts/release_darwin.sh # this script
##                  logs/                     # build/test logs
##                  source/                   # svn checkout of release branch
##                  flavor/opt|debug          # released libraries/headers
## 2. Assumes that this script lives at the location above when
##    executed.

##---------------------------------------------------------------------------##
## Instructions
##---------------------------------------------------------------------------##
## 1. Set modulefiles to be loaded in named environment functions.
## 2. Update variables that control the build:
##    - $ddir
##    - $CONFIG_BASE
## 3. Run this script: ./release_ml &> ../logs/relase_moonlight.log

#----------------------------------------------------------------------#
# Per release settings go here (edits go here)
#----------------------------------------------------------------------#

# Draco install directory name (/usr/projects/draco/draco-NN_NN_NN)
export package=draco
ddir=draco-6_18_0
pdir=$ddir

# CMake options that will be included in the configuration step
export CONFIG_BASE="-DDRACO_VERSION_PATCH=0 -DUSE_CUDA=OFF"

# environment (use draco modules)
# release for each module set
environments="intel15env"
function intel14env()
{
  run "module purge"
  run "module use --append /usr/projects/draco/vendors/Modules"
  run "module load cmake/3.3.2 numdiff/5.2.1 python/2.7.3"
  run "module load compilers/intel/14.0.2 mpi/openmpi-1.6.5-intel_14.0.2"
  run "module load random123 eospac/6.2.4"
  run "module list"
  export MPIEXEC=${MPIRUN}
}
function intel15env()
{
  run "module purge"
  run "module use --append /usr/projects/draco/vendors/Modules"
  run "module load cmake/3.4.0 numdiff/5.8.1"
  run "module load intel/15.0.3 openmpi/1.6.5-intel_15.0.3"
  run "module load random123 eospac/6.2.4"
  run "module list"
  export MPIEXEC=${MPIRUN}
}

# ============================================================================
# ====== Normally, you do not edit anything below this line ==================
# ============================================================================

##---------------------------------------------------------------------------##
## Generic setup (do not edit)
##---------------------------------------------------------------------------##
sdir=`dirname $0`
cdir=`pwd`
cd $sdir
export script_dir=`pwd`
export draco_script_dir=$script_dir
cd $cdir
source $draco_script_dir/common.sh

# sets umask 0002
# sets $install_group, $install_permissions, $build_permissions
establish_permissions

export source_prefix="/usr/projects/$package/$pdir"
scratchdir=`selectscratchdir`

# =============================================================================
# Build types:
# - These must be copied into release_ml.msub because bash arrays cannot
#   be passed to the subshell (bash bug)
# =============================================================================

OPTIMIZE_ON="-DCMAKE_BUILD_TYPE=Release"
OPTIMIZE_OFF="-DCMAKE_BUILD_TYPE=Debug  "

LOGGING_ON="-DDRACO_DIAGNOSTICS=7 -DDRACO_TIMING=1"
LOGGING_OFF="-DDRACO_DIAGNOSTICS=0 -DDRACO_TIMING=0"

# Define the meanings of the various code versions:

VERSIONS=( "debug" "opt" )
OPTIONS=(\
    "$OPTIMIZE_OFF $LOGGING_OFF" \
    "$OPTIMIZE_ON $LOGGING_OFF" \
)

##---------------------------------------------------------------------------##
## Environment review
##---------------------------------------------------------------------------##

verbose=1
if test $verbose == 1; then
  echo
  echo "Build environment summary:"
  echo "=========================="
  echo "script_dir     = $script_dir"
  echo "source_prefix  = $source_prefix"
  echo "log_dir        = $source_prefix/logs"
  echo
  echo "package     = $package"
  echo "versions:"
  for (( i=0 ; i < ${#VERSIONS[@]} ; ++i )); do
    echo -e "   ${VERSIONS[$i]}, \t options = ${OPTIONS[$i]}"
  done
  echo
fi

##---------------------------------------------------------------------------##
## Execute the build, test and install
##---------------------------------------------------------------------------##

jobids=""
for env in $environments; do

  # Run the bash function defined above to load appropriate module
  # environment.
  echo -e "\nEstablish environment $env"
  echo "======================================="
  $env

  buildflavor=`flavor`
  # e.g.: buildflavor=moonlight-openmpi-1.6.5-intel-15.0.3

  export install_prefix="$source_prefix/$buildflavor"
  export build_prefix="$scratchdir/$USER/$pdir/$buildflavor"
  export draco_prefix="/usr/projects/draco/$ddir/$buildflavor"

  for (( i=0 ; i < ${#VERSIONS[@]} ; ++i )); do

    export version=${VERSIONS[$i]}
    export options=${OPTIONS[$i]}

    export CONFIG_EXTRA="$CONFIG_BASE"

    # export dry_run=1
    # https://darwin.lanl.gov/darwin_hw/report.html
    export steps="config build test"
    cmd="/usr/bin/sbatch -v -N 1 -p sl230s -t 360 \
-o $source_prefix/logs/release-$buildflavor-$version.log \
-e $source_prefix/logs/release-$buildflavor-$version.log \
$script_dir/release_darwin.msub"
    echo -e "\nConfigure, Build and Test $buildflavor-$version version of $package."
    echo "$cmd"
    jobid=`eval ${cmd} &`
    sleep 1m
    # trim extra whitespace from number
    jobid=`echo ${jobid//[^0-9]/}`
    export jobids="$jobid $jobids"

    # export dry_run=0

  done
done

##---------------------------------------------------------------------------##
## Set permissions
##---------------------------------------------------------------------------##

publish_release

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
