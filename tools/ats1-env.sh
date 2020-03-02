#!/bin/bash
#------------------------------------------------------------------------------#
# Cray (ATS-1) Environment setups
#------------------------------------------------------------------------------#

# verbose=true
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
if [[ `df | grep yellnfs2 | grep -c jayenne` -gt 0 ]]; then
  export siblings="trinitite"
else
  export siblings="trinity"
fi

# The following toolchains will be used when releasing code
environments="intel1904env intel1904env-knl intel1802env intel1802env-knl intel1704env intel1704env-knl"

# Extra cmake options
export CONFIG_BASE+=" -DCMAKE_VERBOSE_MAKEFILE=ON"

# SLURM
avail_queues=`sacctmgr -np list assoc user=$LOGNAME | sed -e 's/.*|\(.*dev.*\|.*access.*\)|.*/\1/' | sed -e 's/|.*//'`
case $avail_queues in
  *access*) access_queue="-A access --qos=access" ;;
  *dev*) access_queue="--qos=dev" ;;
esac
export access_queue

#------------------------------------------------------------------------------#
# Specify environments (modules)
#------------------------------------------------------------------------------#

if ! [[ $ddir ]] ;then
  echo "FATAL ERROR: Expected ddir to be set in the environment. (ats1-env.sh)"
  exit 1
fi

case $ddir in

  #------------------------------------------------------------------------------#
  draco-7_2* | draco-7_3* | draco-7_4*)
    function intel1904env()
    {
      unset partition
      unset jobnameext

      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi quo"
      run "module unload superlu-dist metis parmetis"
      run "module unload csk lapack"
      run "module unload PrgEnv-intel PrgEnv-pgi PrgEnv-cray PrgEnv-gnu"
      run "module unload lapack "
      run "module unload intel gcc"
      run "module unload papi perftools"
      run "module load PrgEnv-intel"
      run "module unload intel"
      run "module unload xt-libsci xt-totalview"
      run "module unload cray-hugepages2M"
      run "module load intel/19.0.4"
      run "module load cmake/3.14.6 numdiff git"
      run "module load gsl random123 eospac/6.4.0 ndi python/3.6-anaconda-5.0.1"
      run "module load trilinos/12.10.1 metis parmetis/4.0.3 superlu-dist quo"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load csk"
      run "module list"
      CC=`which cc`
      CXX=`which CC`
      FC=`which ftn`
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }
    function intel1802env()
    {
      unset partition
      unset jobnameext

      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib"
      run "module load friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi quo"
      run "module unload superlu-dist metis parmetis"
      run "module unload csk lapack"
      run "module unload PrgEnv-intel PrgEnv-pgi PrgEnv-cray PrgEnv-gnu"
      run "module unload lapack "
      run "module unload intel gcc"
      run "module unload papi perftools"
      run "module load PrgEnv-intel"
      run "module unload intel"
      run "module unload xt-libsci xt-totalview"
      run "module unload cray-hugepages2M"
      run "module load intel/18.0.2"
      run "module load cmake/3.14.6 numdiff git"
      run "module load gsl random123 eospac/6.4.0 ndi python/3.6-anaconda-5.0.1"
      run "module load trilinos/12.10.1 metis parmetis/4.0.3 superlu-dist quo"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load csk"
      run "module list"
      CC=`which cc`
      CXX=`which CC`
      FC=`which ftn`
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }
    function intel1904env-knl()
    {
      export partition="-p knl"
      export jobnameext="-knl"

      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi quo"
      run "module unload superlu-dist metis parmetis"
      run "module unload csk lapack"
      run "module unload PrgEnv-intel PrgEnv-pgi PrgEnv-cray PrgEnv-gnu"
      run "module unload intel gcc"
      run "module unload papi perftools"
      run "module load PrgEnv-intel"
      run "module unload intel"
      run "module unload xt-libsci xt-totalview"
      run "module unload cray-hugepages2M"
      run "module load intel/19.0.4"
      run "module load cmake/3.14.6 numdiff git"
      run "module load gsl random123 eospac/6.4.0 ndi python/3.6-anaconda-5.0.1"
      run "module load trilinos/12.10.1 metis parmetis/4.0.3 superlu-dist quo"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load csk"
      run "module swap craype-haswell craype-mic-knl"
      run "module list"
      run "module list"
      CC=`which cc`
      CXX=`which CC`
      FC=`which ftn`
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }

    function intel1802env-knl()
    {
      export partition="-p knl"
      export jobnameext="-knl"

      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi quo"
      run "module unload superlu-dist metis parmetis"
      run "module unload csk lapack"
      run "module unload PrgEnv-intel PrgEnv-pgi PrgEnv-cray PrgEnv-gnu"
      run "module unload intel gcc"
      run "module unload papi perftools"
      run "module load PrgEnv-intel"
      run "module unload intel"
      run "module unload xt-libsci xt-totalview"
      run "module unload cray-hugepages2M"
      run "module load intel/18.0.2"
      run "module load cmake/3.14.6 numdiff git"
      run "module load gsl random123 eospac/6.4.0 ndi python/3.6-anaconda-5.0.1"
      run "module load trilinos/12.10.1 metis parmetis/4.0.3 superlu-dist quo"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load csk"
      run "module swap craype-haswell craype-mic-knl"
      run "module list"
      run "module list"
      CC=`which cc`
      CXX=`which CC`
      FC=`which ftn`
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }
    function intel1704env()
    {
      unset partition
      unset jobnameext

      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi quo"
      run "module unload superlu-dist metis parmetis"
      run "module unload csk lapack"
      run "module unload PrgEnv-intel PrgEnv-pgi PrgEnv-cray PrgEnv-gnu"
      run "module unload lapack "
      run "module unload intel gcc"
      run "module unload papi perftools"
      run "module load PrgEnv-intel"
      run "module unload intel"
      run "module unload xt-libsci xt-totalview"
      run "module load intel/17.0.4"
      run "module load cmake/3.14.6 numdiff git"
      run "module load gsl random123 eospac/6.4.0 ndi python/3.6-anaconda-5.0.1"
      run "module load trilinos/12.10.1 metis parmetis/4.0.3 superlu-dist quo"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load csk"
      run "module list"
      CC=`which cc`
      CXX=`which CC`
      FC=`which ftn`
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }
    function intel1704env-knl()
    {
      export partition="-p knl"
      export jobnameext="-knl"

      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi quo"
      run "module unload superlu-dist metis parmetis"
      run "module unload csk lapack"
      run "module unload PrgEnv-intel PrgEnv-pgi PrgEnv-cray PrgEnv-gnu"
      run "module unload intel gcc"
      run "module unload papi perftools"
      run "module load PrgEnv-intel"
      run "module unload intel"
      run "module unload xt-libsci xt-totalview"
      run "module load intel/17.0.4"
      run "module load cmake/3.14.6 numdiff git"
      run "module load gsl random123 eospac/6.4.0 ndi python/3.6-anaconda-5.0.1"
      run "module load trilinos/12.10.1 metis parmetis/4.0.3 superlu-dist quo"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load csk"
      run "module swap craype-haswell craype-mic-knl"
      run "module list"
      run "module list"
      CC=`which cc`
      CXX=`which CC`
      FC=`which ftn`
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }
    ;;

  #------------------------------------------------------------------------------#
  draco-7_0* | draco-7_1*)
    function intel1802env()
    {
      unset partition
      unset jobnameext

      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi quo"
      run "module unload superlu-dist metis parmetis"
      run "module unload csk lapack"
      run "module unload PrgEnv-intel PrgEnv-pgi PrgEnv-cray PrgEnv-gnu"
      run "module unload lapack "
      run "module unload intel gcc"
      run "module unload papi perftools"
      run "module load PrgEnv-intel"
      run "module unload intel"
      run "module unload xt-libsci xt-totalview"
      run "module unload cray-hugepages2M"
      run "module load intel/18.0.2"
      run "module load cmake/3.12.1 numdiff git"
      run "module load gsl random123 eospac/6.3.0 ndi"
      run "module load trilinos/12.10.1 metis parmetis/4.0.3 superlu-dist"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load csk"
      run "module list"
      CC=`which cc`
      CXX=`which CC`
      FC=`which ftn`
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }
    function intel1802env-knl()
    {
      export partition="-p knl"
      export jobnameext="-knl"
      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi quo"
      run "module unload superlu-dist metis parmetis"
      run "module unload csk lapack"
      run "module unload PrgEnv-intel PrgEnv-pgi PrgEnv-cray PrgEnv-gnu"
      run "module unload intel gcc"
      run "module unload papi perftools"
      run "module load PrgEnv-intel"
      run "module unload intel"
      run "module unload xt-libsci xt-totalview"
      run "module unload cray-hugepages2M"
      run "module load intel/18.0.2"
      run "module load cmake/3.12.1 numdiff git"
      run "module load gsl random123 eospac/6.3.0 ndi"
      run "module load trilinos/12.10.1 metis parmetis/4.0.3 superlu-dist"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load csk"
      run "module swap craype-haswell craype-mic-knl"
      run "module list"
      run "module list"
      CC=`which cc`
      CXX=`which CC`
      FC=`which ftn`
      export TARGET=knl
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }
    function intel1704env()
    {
      unset partition
      unset jobnameext

      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi quo"
      run "module unload superlu-dist metis parmetis"
      run "module unload csk lapack"
      run "module unload PrgEnv-intel PrgEnv-pgi PrgEnv-cray PrgEnv-gnu"
      run "module unload lapack "
      run "module unload intel gcc"
      run "module unload papi perftools"
      run "module load PrgEnv-intel"
      run "module unload intel"
      run "module unload xt-libsci xt-totalview"
      run "module load intel/17.0.4"
      run "module load cmake/3.12.1 numdiff git"
      run "module load gsl random123 eospac/6.3.0 ndi"
      run "module load trilinos/12.10.1 metis parmetis/4.0.3 superlu-dist"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load csk"
      run "module list"
      CC=`which cc`
      CXX=`which CC`
      FC=`which ftn`
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }
    function intel1704env-knl()
    {
      export partition="-p knl"
      export jobnameext="-knl"
      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi quo"
      run "module unload superlu-dist metis parmetis"
      run "module unload csk lapack"
      run "module unload PrgEnv-intel PrgEnv-pgi PrgEnv-cray PrgEnv-gnu"
      run "module unload intel gcc"
      run "module unload papi perftools"
      run "module load PrgEnv-intel"
      run "module unload intel"
      run "module unload xt-libsci xt-totalview"
      run "module load intel/17.0.4"
      run "module load cmake/3.12.1 numdiff git"
      run "module load gsl random123 eospac/6.3.0 ndi"
      run "module load trilinos/12.10.1 metis parmetis/4.0.3 superlu-dist"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load csk"
      run "module swap craype-haswell craype-mic-knl"
      run "module list"
      run "module list"
      CC=`which cc`
      CXX=`which CC`
      FC=`which ftn`
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }
    ;;

  *)
    die "ats1-env.sh:: did not set any build environments, ddir = $ddir."
    ;;


  #------------------------------------------------------------------------------#
esac

#------------------------------------------------------------------------------#
# Sanity check
#------------------------------------------------------------------------------#

for env in $environments; do
  if [[ `fn_exists $env` -gt 0 ]]; then
    if [[ $verbose ]]; then echo "export -f $env"; fi
    # exporting these environment is required on CTS-1, but breaks the module
    # system on trinitite.  Ugh.
    # export -f $env
  else
    die "Requested environment $env is not defined."
  fi
done

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
