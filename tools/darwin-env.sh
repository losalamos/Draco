#!/bin/bash
#-----------------------------------------------------------------------------#
# Darwin Environment setups
#-----------------------------------------------------------------------------#

if [[ `uname -n` =~ "darwin-fe" ]]; then
  echo "FATAL ERROR: This script must be run from a back-end node"
  exit 1
fi

if ! [[ $ddir ]] ;then
  echo "FATAL ERROR: Expected ddir to be set in the environment. (darwin-env.sh)"
  exit 1
fi

export VENDOR_DIR=/projects/draco/vendors
export DRACO_ARCH=`/usr/projects/draco/vendors/bin/target_arch`

# The following toolchains will be used when releasing code:
# case ${SLURM_JOB_PARTITION:-unknown} in
#   power9*)   environments="p9gcc730env" ;;
#   volta*x86) environments="x86gcc730env" ;;
#   arm*)      environments="armgcc820env" ;;
# esac

# Extra cmake options
export CONFIG_BASE+=" -DCMAKE_VERBOSE_MAKEFILE=ON"

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
