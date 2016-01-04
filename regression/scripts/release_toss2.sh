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
# Per release settings go here (edits go here)
#----------------------------------------------------------------------#

# Draco install directory name (/usr/projects/draco/draco-NN_NN_NN)
export package=draco
ddir=draco-6_19_0
pdir=$ddir

# CMake options that will be included in the configuration step
export CONFIG_BASE="-DDRACO_VERSION_PATCH=0"

# environment (use draco modules)
# release for each module set
environments="intel15env"
function intel14env()
{
  run "module use /usr/projects/draco/vendors/Modules/hpc"
  run "module purge"
  run "module load friendly-testing user_contrib"
  run "module load cmake/3.4.0 svn numdiff"
  run "module load intel/14.0.4 openmpi/1.6.5"
  run "module load random123 eospac/6.2.4"
  run "module list"
}
function intel15env()
{
  run "module use /usr/projects/draco/vendors/Modules/hpc"
  run "module purge"
  run "module load friendly-testing user_contrib"
  run "module load cmake/3.4.0 svn numdiff"
  run "module load intel/15.0.5 openmpi/1.6.5"
  run "module load random123 eospac/6.2.4"
  run "module list"
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
ppn=`showstats -n | grep 10 | head -n 1 | awk '{print $3}'`
#build_pe=`npes_build`
#test_pe=`npes_test`

avail_queues=`mdiag -u $LOGNAME | grep ALIST | sed -e 's/.*ALIST=//' | sed -e 's/,/ /g'`
case $avail_queues in
  *access*) access_queue="-A access" ;;
esac

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
  echo "scratchdir     = /$scratchdir/$USER"
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
  export build_prefix="/$scratchdir/$USER/$pdir/$buildflavor"
  export draco_prefix="/usr/projects/draco/$ddir/$buildflavor"

  for (( i=0 ; i < ${#VERSIONS[@]} ; ++i )); do

    export version=${VERSIONS[$i]}
    export options=${OPTIONS[$i]}

    export CONFIG_EXTRA="$CONFIG_BASE"

    # export dry_run=1
    export steps="config build test"
    cmd="msub -V $access_queue -l walltime=01:00:00 -l nodes=1:ppn=${ppn} -j oe \
-o $source_prefix/logs/release-$buildflavor-$version.log $script_dir/release_toss2.msub"
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

# The Pack_Build_gnu EAP target needs this symlink on moonlight
if test `machineName` == moonlight; then
  run "cd $source_prefix"
  gccflavor=`echo $flavor | sed -e s%$LMPI-$LMPIVER%gcc-4.9.2%`
  run "ln -s $flavor $gccflavor"
fi

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
