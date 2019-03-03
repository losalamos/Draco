#!/bin/bash
#------------------------------------------------------------------------------#
# Cray (ATS-1) Environment setups
#------------------------------------------------------------------------------#

export VENDOR_DIR=/usr/projects/draco/vendors
if [[ -d $ParMETIS_ROOT_DIR ]]; then
  echo "ERROR: This script should be run from a clean environment."
  echo "       Try running 'rmdracoenv'."
  exit 1
fi

case $ddir in

  #------------------------------------------------------------------------------#
  draco-7_0_0 | draco-7_1_0)
    function intel1802env()
    {
      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi"
      run "module unload superlu-dist metis parmetis"
      run "module unload csk lapack"
      run "module unload PrgEnv-intel PrgEnv-pgi PrgEnv-cray PrgEnv-gnu"
      run "module unload lapack "
      run "module unload intel gcc"
      run "module unload papi perftools"
      run "module load PrgEnv-intel"
      run "module unload intel"
      run "module unload xt-libsci xt-totalview"
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
      export CRAYPE_LINK_TYPE=dynamic
      export OMP_NUM_THREADS=16
      export TARGET=haswell
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }
    function intel1802env-knl()
    {
      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi"
      run "module unload superlu-dist metis parmetis"
      run "module unload csk lapack"
      run "module unload PrgEnv-intel PrgEnv-pgi PrgEnv-cray PrgEnv-gnu"
      run "module unload intel gcc"
      run "module unload papi perftools"
      run "module load PrgEnv-intel"
      run "module unload intel"
      run "module unload xt-libsci xt-totalview"
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
      export CRAYPE_LINK_TYPE=dynamic
      export OMP_NUM_THREADS=17
      export TARGET=knl
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }
    function intel1704env()
    {
      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi"
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
      export CRAYPE_LINK_TYPE=dynamic
      export OMP_NUM_THREADS=16
      export TARGET=haswell
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }
    function intel1704env-knl()
    {
      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi"
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
      export CRAYPE_LINK_TYPE=dynamic
      export OMP_NUM_THREADS=17
      export TARGET=knl
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }
    ;;

    #------------------------------------------------------------------------------#
    draco-6_25_0 )

    function intel18env()
    {
      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi"
      run "module unload superlu-dist metis parmetis"
      run "module unload csk lapack"
      run "module unload PrgEnv-intel PrgEnv-pgi PrgEnv-cray PrgEnv-gnu"
      run "module unload lapack "
      run "module unload intel gcc"
      run "module unload papi perftools"
      run "module load PrgEnv-intel"
      run "module unload intel"
      run "module unload xt-libsci xt-totalview"
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
      export CRAYPE_LINK_TYPE=dynamic
      export OMP_NUM_THREADS=16
      export TARGET=haswell
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }

    function intel18env-knl()
    {
      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi"
      run "module unload superlu-dist metis parmetis"
      run "module unload csk lapack"
      run "module unload PrgEnv-intel PrgEnv-pgi PrgEnv-cray PrgEnv-gnu"
      run "module unload intel gcc"
      run "module unload papi perftools"
      run "module load PrgEnv-intel"
      run "module unload intel"
      run "module unload xt-libsci xt-totalview"
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
      export CRAYPE_LINK_TYPE=dynamic
      export OMP_NUM_THREADS=17
      export TARGET=knl
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }

    function intel17env()
    {
      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi"
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
      export CRAYPE_LINK_TYPE=dynamic
      export OMP_NUM_THREADS=16
      export TARGET=haswell
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }

    function intel17env-knl()
    {
      if [[ ${CRAY_CPU_TARGET} == mic-knl ]]; then
        run "module swap craype-mic-knl craype-haswell"
      fi
      run "module load user_contrib friendly-testing"
      run "module unload cmake numdiff git"
      run "module unload gsl random123 eospac"
      run "module unload trilinos ndi"
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
      export CRAYPE_LINK_TYPE=dynamic
      export OMP_NUM_THREADS=17
      export TARGET=knl
      export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    }
    ;;


  #------------------------------------------------------------------------------#
esac


##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
