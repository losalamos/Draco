#!/bin/bash
#------------------------------------------------------------------------------#
# Cray (ARM) Environment setups (capulin, thunder)
#------------------------------------------------------------------------------#

export KMP_AFFINITY=disabled
export VENDOR_DIR=/usr/projects/draco/vendors
export CRAYPE_LINK_TYPE=dynamic

# Sanity Check (Cray machines have very fragile module systems!)
if [[ -d $ParMETIS_ROOT_DIR ]]; then
  echo "ERROR: This script should be run from a clean environment."
  echo "       Try running 'rmdracoenv'."
  exit 1
fi

# symlinks will be generated for each machine that point to the correct
# installation directory.
if [[ `df | grep yellow | grep -c lustre` -gt 0 ]]; then
  export siblings="capulin"
else
  export siblings="thunder"
fi

# The following toolchains will be used when releasing code
environments="gcc830 cce90180"

# Extra cmake options
export CONFIG_BASE+=" -DCMAKE_VERBOSE_MAKEFILE=ON"

#------------------------------------------------------------------------------#
# Specify environments (modules)
#------------------------------------------------------------------------------#

if ! [[ $ddir ]] ;then
  echo "FATAL ERROR: Expected ddir to be set in the environment. (cray-env.sh)"
  exit 1
fi

case $ddir in

  #------------------------------------------------------------------------------#
  draco-7_3* | draco-7_4*)
    function gcc830()
    {
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi"
      run "module unload superlu-dist metis parmetis"
      run "module unload csk lapack"
      run "module unload PrgEnv-allinea PrgEnv-cray PrgEnv-gnu"
      run "module unload lapack "
      run "module unload cce gcc"
      run "module unload papi perftools"
      run "module load PrgEnv-gnu"
      run "module load cmake numdiff/5.9.0-gcc-8.3.0 gsl/2.5-gcc-8.3.0"
      run "module load openblas/0.3.6-gcc-8.3.0 metis/5.1.0-gcc-8.3.0"
      run "module load eospac/6.4.0-gcc-8.3.0 random123/1.09-gcc-8.3.0"
      run "module load parmetis/4.0.3-gcc-8.3.0-mpt-7.7.8"
      run "module load superlu-dist/5.4.0-gcc-8.3.0-mpt-7.7.8-openblas"
      run "module load trilinos/12.14.1-gcc-8.3.0-mpt-7.7.8-openblas"
      run "module load cray-python/3.6.5.6 qt"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles-capulin"
      run "module load csk/0.4.2-gcc-8.3.0"
      run "module list"
      export CC=`which cc`
      export CXX=`which CC`
      export FC=`which ftn`
#      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }

    function cce90180()
    {
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi"
      run "module unload superlu-dist metis parmetis"
      run "module unload csk lapack"
      run "module unload PrgEnv-allinea PrgEnv-cray PrgEnv-gnu"
      run "module unload lapack "
      run "module unload cce gcc"
      run "module unload papi perftools"
      run "module load PrgEnv-cray"
      run "module load cmake numdiff/5.9.0-cce-9.0.1.80 gsl/2.5-cce-9.0.1.80"
      run "module load metis/5.1.0-cce-9.0.1.80"
      run "module load eospac/6.4.0-cce-9.0.1.80 random123/1.09-cce-9.0.1.80"
      run "module load parmetis/4.0.3-cce-9.0.1.80-mpt-7.7.8"
      # run "module load superlu-dist/5.4.0-cce-9.0.1.80-mpt-7.7.8-openblas"
      run "module load cray-python/3.6.5.6 qt"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles-capulin"
      run "module load csk/0.4.2-cce-9.0.1.80"
      run "module list"
      export CC=`which cc`
      export CXX=`which CC`
      export FC=`which ftn`
#      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }

    ;;

  *)
    die "cray-env.sh:: did not set any build environments, ddir = $ddir."
    ;;

  #------------------------------------------------------------------------------#
esac

#------------------------------------------------------------------------------#
# Sanity check
#------------------------------------------------------------------------------#

for env in $environments; do
  if [[ `fn_exists $env` -gt 0 ]]; then
    if [[ $verbose ]]; then echo "export -f $env"; fi
    export -f $env
  else
    die "Requested environment $env is not defined."
  fi
done


##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
