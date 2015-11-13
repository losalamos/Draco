#!/bin/bash -l

##---------------------------------------------------------------------------##
## Assumptions:
##---------------------------------------------------------------------------##
## 1. Directory layout:
##    /usr/projects/draco/draco-NN_NN_NN/
##                  scripts/release_toss2.sh # this script
##                  logs/                    # build/test logs
##                  source/draco-NN_NN_NN    # svn checkout of release branch
##                  flavor/opt|debug         # released libraries/headers
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
# Per release settings go here:
#----------------------------------------------------------------------#

# Draco install directory name (/usr/projects/draco/draco-NN_NN_NN)
export package=draco
ddir=draco-6_18_0
pdir=$ddir

# CMake options that will be included in the configuration step
export CONFIG_BASE="-DDRACO_VERSION_PATCH=0"

# environment (use draco modules)
# release for each module set
target="`uname -n | sed -e s/[.].*//`"
case $target in
  c[it]-fe[0-9] | c[it]-login[0-9] | c[it]-vizlogin[0-9])
    environments="intel14env intel15env" ;;
  tt-fey* | tt-login*)
    environments="intel15env" ;;
esac
function intel15env()
{
run "module load friendly-testing user_contrib"
run "module unload ndi ParMetis SuperLU_DIST trilinos"
run "module unload lapack gsl intel"
run "module unload cmake numdiff svn"
run "module unload PrgEnv-intel PrgEnv-pgi PrgEnv-cray PrgEnv-gnu"
run "module unload papi perftools"
run "module load PrgEnv-intel"
run "module unload xt-libsci xt-totalview"
case $target in
  c[it]-fe[0-9] | c[it]-login[0-9] | c[it]-vizlogin[0-9])
    run "module swap intel intel/15.0.3"
    run "module load gsl/1.15"
    ;;
  tt-fey* | tt-login*)
    run "module swap intel intel/15.0.3.187"
    run "module load gsl/1.16"
    ;;
esac
run "module load cmake/3.3.2 numdiff svn"
run "module load trilinos SuperLU_DIST"
run "module load ParMetis ndi random123 eospac/6.2.4"
run "module list"
CC=`which cc`
CXX=`which CC`
FC=`which ftn`
export OMP_NUM_THREADS=8
}

function intel14env()
{
run "module load friendly-testing user_contrib"
run "module unload ndi ParMetis SuperLU_DIST trilinos"
run "module unload lapack gsl intel"
run "module unload cmake numdiff svn"
run "module unload PrgEnv-intel PrgEnv-pgi"
run "module unload papi perftools"
run "module load PrgEnv-intel"
run "module unload xt-libsci xt-totalview"
run "module swap intel intel/14.0.4.211"
run "module load gsl/1.15"
run "module load cmake/3.3.2 numdiff svn"
run "module load trilinos SuperLU_DIST"
run "module load ParMetis ndi random123 eospac/6.2.4"
run "module list"
CC=`which cc`
CXX=`which CC`
FC=`which ftn`
export OMP_NUM_THREADS=8
}

# ============================================================================
# ====== Normally, you do not edit anything below this line ==================
# ============================================================================

##---------------------------------------------------------------------------##
## Generic setup
##---------------------------------------------------------------------------##
sdir=`dirname $0`
cdir=`pwd`
cd $sdir
export script_dir=`pwd`
export draco_script_dir=$script_dir
cd $cdir
source $script_dir/common.sh

# sets umask 0002
# sets $install_group, $install_permissions, $build_permissions
establish_permissions

export source_prefix="/usr/projects/$package/$pdir"
scratchdir=`selectscratchdir`
ppn=`lookupppn`

# =============================================================================
# Build types:
# - These must be copied into release_ml.msub because bash arrays cannot
#   be passed to the subshell (bash bug)
# =============================================================================

OPTIMIZE_ON="-DCMAKE_BUILD_TYPE=Release -DDRACO_LIBRARY_TYPE=STATIC"
OPTIMIZE_OFF="-DCMAKE_BUILD_TYPE=Debug  -DDRACO_LIBRARY_TYPE=STATIC"
OPTIMIZE_RWDI="-DCMAKE_BUILD_TYPE=RelWithDebInfo -DDRACO_LIBRARY_TYPE=SHARED"

LOGGING_ON="-DDRACO_DIAGNOSTICS=7 -DDRACO_TIMING=1"
LOGGING_OFF="-DDRACO_DIAGNOSTICS=0 -DDRACO_TIMING=0"

# Define the meanings of the various code versions:

VERSIONS=( "debug" "opt" "rwdi" )
OPTIONS=(\
    "$OPTIMIZE_OFF  $LOGGING_OFF" \
    "$OPTIMIZE_ON   $LOGGING_OFF" \
    "$OPTIMIZE_RWDI $LOGGING_OFF" \
)

##---------------------------------------------------------------------------##
## Environment review
##---------------------------------------------------------------------------##

verbose=1
if test $verbose == 1; then
  echo
  echo "Build environment summary:"
  echo "=========================="
  echo "script_dir       = $script_dir"
  echo "draco_script_dir = $script_dir"
  echo "source_prefix    = $source_prefix"
  echo "log_dir          = $source_prefix/logs"
  echo
  echo "package          = $package"
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
  export build_prefix="/$scratchdir/$USER/$pdir/$buildflavor"
  export draco_prefix="/usr/projects/draco/$ddir/$buildflavor"

  for (( i=0 ; i < ${#VERSIONS[@]} ; ++i )); do

    export version=${VERSIONS[$i]}
    export options=${OPTIONS[$i]}

    export CONFIG_EXTRA="$CONFIG_BASE"

    # export dry_run=1
    # config and build on front-end
    echo -e "\nConfigure and build $package for $buildflavor-$version."
    export steps="config build"
    run "$draco_script_dir/release_cray.msub &> $source_prefix/logs/release-$buildflavor-$version-cb.log"

    # Run the tests on the back-end.
    export steps="test"
    cmd="msub -V $access_queue -l walltime=08:00:00 -l nodes=1:ppn=${ppn} -j oe \
-o $source_prefix/logs/release-$buildflavor-$version-t.log $draco_script_dir/release_cray.msub"
    echo -e "\nTest $package for $buildflavor-$version."
    echo "$cmd"
    jobid=`eval ${cmd}`
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
