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

# Avoid 'unbound variable' errors on Sequoia
set +u

# CAUTION: Cannot have too many environment variables set!  When
# running this script use a bare-bones environment and unset things
# like LS_COLORS and MANPATH.

# Draco install directory name (/usr/projects/draco/draco-NN_NN_NN)
export package=draco
ddir=draco-6_20_1
pdir=$ddir

# environment (use draco modules)
# release for each module set
environments="gcc484"
function gcc484()
{
  export VENDOR_DIR=/usr/gapps/jayenne/vendors
  export DK_NODE=$DK_NODE:/$VENDOR_DIR/Modules/sq
  export OMP_NUM_THREADS=4
  use gcc484 python-2.7.3
  use cmake361 gsl numdiff random123
  use
}
function xlc12()
{
  export VENDOR_DIR=/usr/gapps/jayenne/vendors
  export DK_NODE=$DK_NODE:/$VENDOR_DIR/Modules/sq
  export OMP_NUM_THREADS=4
  use xlc12
  use cmake340 gsl numdiff random123
  use
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

# CMake options that will be included in the configuration step
export CONFIG_BASE="-DDRACO_LIBRARY_TYPE=STATIC -DDRACO_VERSION_PATCH=`echo $ddir | sed -e 's/.*_//'`"

# sets umask 0002
# sets $install_group, $install_permissions, $build_permissions
establish_permissions

export source_prefix="/usr/gapps/jayenne/$pdir"
scratchdir=`selectscratchdir`
build_pe=`npes_build`
test_pe=`npes_test`

# Sequoia limits the environment size to 8192, so clean up to keep the
# environment below this limit.
unset LS_COLORS
unset MANPATH
unset NLSPATH
unset EXINIT
unset CLASSPATH
unset QTDIR
unset QTLIB
unset QTINC
unset DRACO_AUTO_CLANG_FORMAT
unset EDITOR
unset SSH_ASKPASS
unset CVS_RSH
unset CEI_HOME

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
  echo "scratchdir     = $scratchdir/$USER"
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
  export draco_prefix="/usr/gapps/jayenne/$ddir/$buildflavor"

  for (( i=0 ; i < ${#VERSIONS[@]} ; ++i )); do

    export version=${VERSIONS[$i]}
    export options=${OPTIONS[$i]}

    export CONFIG_EXTRA="$CONFIG_BASE"

    # export dry_run=1

    # Config and build on the front end
    echo -e "\nConfigure and Build $buildflavor-$version version of $package."
    export steps="config build"
    run "$draco_script_dir/release_bgq.msub &> $source_prefix/logs/release-$buildflavor-$version-cb.log"

    # Run the tests on the back-end.
    echo -e "\nTest $buildflavor-$version version of $package."
    export steps="test"
    cmd="sbatch -v -N 64 -p pdebug -t 6:00:00 \
-o $source_prefix/logs/release-$buildflavor-$version-t.log \
-e $source_prefix/logs/release-$buildflavor-$version-t.log \
$script_dir/release_bgq.msub"
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
