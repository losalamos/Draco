#!/bin/bash
#------------------------------------------------------------------------------#
# CTS-1 Environment setups
#------------------------------------------------------------------------------#

case $ddir in

  #------------------------------------------------------------------------------#
  draco-6_25_0 | draco-7_0_0 | draco-7_1_0)
    function intel1802env()
    {
      export VENDOR_DIR=/usr/projects/draco/vendors
      run "module purge"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load friendly-testing user_contrib"
      run "module load cmake git numdiff python/3.6-anaconda-5.0.1"
      run "module load intel/18.0.2 openmpi/2.1.2"
      run "unset MPI_ROOT"
      run "module load random123 eospac/6.3.0 gsl"
      run "module load mkl metis ndi csk qt"
      run "module load parmetis superlu-dist trilinos"
      run "module list"
    }
    function intel1704env()
    {
      export VENDOR_DIR=/usr/projects/draco/vendors
      run "module purge"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load friendly-testing user_contrib"
      run "module load cmake git numdiff python/3.6-anaconda-5.0.1"
      run "module load intel/17.0.4 openmpi/2.1.2"
      run "unset MPI_ROOT"
      run "module load random123 eospac/6.3.0 gsl"
      run "module load mkl metis ndi csk qt"
      run "module load parmetis superlu-dist trilinos"
      run "module list"
    }
    function gcc640env()
    {
      export VENDOR_DIR=/usr/projects/draco/vendors
      run "module purge"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load friendly-testing user_contrib"
      run "module load cmake git numdiff python/3.6-anaconda-5.0.1"
      run "module load gcc/6.4.0 openmpi/2.1.2"
      run "unset MPI_ROOT"
      run "module load random123 eospac/6.3.0 gsl"
      run "module load mkl metis ndi qt"
      run "module load parmetis superlu-dist trilinos"
      run "module list"
    }
    ;;

#------------------------------------------------------------------------------#
  draco-6_23_0 )

    function intel1704env()
    {
      run "module purge"
      run "module load friendly-testing user_contrib"
      run "module load cmake git numdiff"
      run "module load intel/17.0.4 openmpi/2.1.2"
      run "module load random123 eospac/6.2.4 gsl"
      run "module load mkl metis ndi csk"
      run "module load parmetis superlu-dist trilinos"
      run "module list"
    }

    function intel1701env()
    {
      run "module purge"
      run "module load friendly-testing user_contrib"
      run "module load cmake git numdiff"
      run "module load intel/17.0.1 openmpi/1.10.5"
      run "module load random123 eospac/6.2.4 gsl"
      run "module load mkl metis ndi csk"
      run "module load parmetis superlu-dist trilinos"
      run "module list"
    }

    function gcc640env()
    {
      run "module purge"
      run "module load friendly-testing user_contrib"
      run "module load cmake git numdiff"
      run "module load gcc/6.4.0 openmpi/2.1.2"
      run "module load random123 eospac/6.2.4 gsl"
      run "module load mkl metis ndi"
      run "module load parmetis superlu-dist trilinos"
      run "module list"
    }
    ;;


  #------------------------------------------------------------------------------#
esac

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
