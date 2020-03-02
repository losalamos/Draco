#!/bin/bash -l

##---------------------------------------------------------------------------##
## Assumptions:
##---------------------------------------------------------------------------##
## 1. Directory layout:
##    /usr/projects/draco/draco-NN_NN_NN/
##                  source/tools/release.sh  # this script
##                  logs/                    # build/test logs
##                  source/                  # git checkout of release branch
##                  flavor/opt|debug|rwdi    # released libraries/headers
## 2. Assumes that this script lives at the location above when executed.
##
## This script pulls data from:
## - <family>-env.sh
## - draco-cmake-opts.sh
##
## This script offloads work to:
## - <family>-release.sh
## - release.msub (via SLURM)

##---------------------------------------------------------------------------##
## Instructions
##---------------------------------------------------------------------------##
## 1. Set modulefiles to be loaded in named environment functions (*-env.sh)
##    - environments (modules, extra variables)
##    - machine list for this family
##    - CONFIG_BASE (extra CMake options) or edit draco-cmake-opts.cmake
## 2. Run this script:
##    cd /usr/projects/draco/draco-NN_NN_NN
##    source/tools/release.sh &> logs/release-snow.log

##---------------------------------------------------------------------------##
## Environment
##---------------------------------------------------------------------------##

# switch to group 'ccsrad' or dacodes
install_group=$USER
build_group=$USER
if [[   `groups | grep -c dacodes` -gt 0 ]]; then install_group=dacodes;
elif [[ `groups | grep -c ccsrad`  -gt 0 ]]; then install_group=ccsrad; fi
if [[ $(id -gn) != $install_group ]]; then exec sg $install_group "$0 $*"; fi
umask 0007       # initially no world or group access
set -m           # Enable job control
shopt -s extglob # Allow variable as case condition
install_permissions="g+rwX,o-rwX"
build_permisssions="g-rwX,o-rwX"

#----------------------------------------------------------------------#
# Per release settings go here (edits go here)
#----------------------------------------------------------------------#

# NOTE: Most settings are specified in the *-env.sh files

# ============================================================================
# ====== Normally, you do not edit anything below this line ==================
# ============================================================================

##---------------------------------------------------------------------------##
## Generic setup (do not edit)
##---------------------------------------------------------------------------##

draco_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
draco_script_dir=`echo $draco_script_dir | awk '{ print $1 }'`
export draco_script_dir=`readlink -f $draco_script_dir`
if [[ -f $draco_script_dir/common.sh ]]; then
  echo "source $draco_script_dir/common.sh"
  source $draco_script_dir/common.sh
else
  echo " "
  echo "FATAL ERROR: Unable to locate Draco's bash functions: "
  echo "   looking for ... source/tools/common.sh"
  echo "   searched draco_script_dir = $draco_script_dir"
  exit 1
fi

export draco_source_prefix=`readlink -f $draco_script_dir/../..`
export ddir=`echo $draco_source_prefix | sed -e 's%.*/%%'`

# If package is set, use those values, otherwise, setup stuff for draco.
if ! [[ $package ]]; then
  export source_prefix=$draco_source_prefix
  export script_dir=$draco_script_dir
  export pdir=$ddir
  export package=`echo $pdir | sed -e 's/-.*//'`
  if [[ -f $draco_script_dir/draco-cmake-opts.sh ]]; then
    echo "source $draco_script_dir/draco-cmake-opts.sh"
    source $draco_script_dir/draco-cmake-opts.sh
  else
    die "Unable to find environment file draco-cmake-opts.sh"
  fi
  # CMake options that will be included in the configuration step
  export CONFIG_BASE+=" -DDRACO_VERSION_PATCH=`echo $ddir | sed -e 's/.*_//'`"
fi

scratchdir=`selectscratchdir`

# hw_threads=`lscpu | grep CPU | head -n 1 | awk '{ print $2 }'`
# hw_threads_per_core=`lscpu | grep Thread | awk '{ print $4 }'`
# ppn=`expr $hw_threads / $hw_threads_per_core`
# build_pe=`npes_build`
# test_pe=`npes_test`

machfam=`machineFamily`
if [[ -f $draco_script_dir/${machfam}-env.sh ]]; then
  echo "source ${draco_script_dir}/${machfam}-env.sh"
  source $draco_script_dir/${machfam}-env.sh
else
  die "Unable to find environment file ${machfam}-env.sh/"
fi

# =============================================================================
# Build types:
# - These must be copied into release.msub because bash arrays cannot be passed
#   to the subshell (bash bug)
# - These variables are set in draco-cmake-opts.sh
# =============================================================================

VERSIONS=( "debug" "opt" "rwdi" )
OPTIONS=(\
    "$OPTIMIZE_OFF $LOGGING_OFF" \
    "$OPTIMIZE_ON $LOGGING_OFF" \
    "$OPTIMIZE_RWDI $LOGGING_OFF" \
)

echo -e "\nBuild environment summary:\n=============================="
echo "script_dir       = $script_dir"
echo "draco_script_dir = $draco_script_dir"
echo "source_prefix    = $source_prefix"
echo "log_dir          = $source_prefix/logs"
echo "scratchdir       = $scratchdir/$USER"
echo -e "\npackage          = $package"
echo "versions:"
for (( i=0 ; i < ${#VERSIONS[@]} ; ++i )); do
  echo -e "   ${VERSIONS[$i]}, \t options = ${OPTIONS[$i]}"; done
echo " "

##---------------------------------------------------------------------------##
## Execute the build, test and install
##---------------------------------------------------------------------------##

jobids=""
echo -e "\nThe following environments will be processed: $environments\n"

for env in $environments; do

  # Run the bash function defined above to load appropriate module
  # environment.
  echo -e "\nEstablish environment $env"
  echo "======================================="
  ${env}

  buildflavor=`flavor`
  # e.g.: buildflavor=snow-openmpi-1.6.5-intel-15.0.3

  export install_prefix="$source_prefix/$buildflavor"
  export build_prefix="$scratchdir/$USER/$pdir/$buildflavor"


  for (( i=0 ; i < ${#VERSIONS[@]} ; ++i )); do

    export rttversion=${VERSIONS[$i]}
    export options=${OPTIONS[$i]}

    echo -e "\nExtra environment setup:\n"

    # callback to append extra data to the cmake options.
    if [[ `fn_exists append_config_base` -gt 0 ]]; then
      extra=`append_config_base`
      echo "  - CONFIG_EXTRA += ${extra}"
      export CONFIG_EXTRA="$extra $CONFIG_BASE"
    else
      export CONFIG_EXTRA="$CONFIG_BASE"
    fi

    source ${draco_script_dir}/${machfam}-release.sh

  done
done

##---------------------------------------------------------------------------##
## Set permissions and set symlinks
##---------------------------------------------------------------------------##

publish_release
wait
relcreatesymlinks

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
