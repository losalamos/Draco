#!/bin/bash
#------------------------------------------------------------------------------#
# ATS-2 Environment setups
#------------------------------------------------------------------------------#

case $ddir in

  #------------------------------------------------------------------------------#
  draco-7_2_0)
    function gcc731env()
    {
      export VENDOR_DIR=/usr/gapps/jayenne/vendors
      run "module purge"
      module use /usr/gapps/jayenne/vendors-ec/spack.20190616/share/spack/lmod/linux-rhel7-ppc64le/Core
      run "module load cuda python/3.7.2 gcc/7.3.1 spectrum-mpi/2019.04.19"
      run "module load cmake/3.14.5 git gsl numdiff random123 metis netlib-lapack"
      run "module load parmetis superlu-dist trilinos csk/0.4.2"
      run "module load eospac/6.4.0"
      # ndi
      run "module list"
      unset MPI_ROOT
      CXX=`which g++`
      CC=`which gcc`
      FC=`which gfortran`
    }
    ;;

  #------------------------------------------------------------------------------#
  draco-7_1_0)
    function gcc731env()
    {
      export VENDOR_DIR=/usr/projects/draco/vendors
      run "module purge"
      run "module use /usr/gapps/user_contrib/spack.20190314/share/spack/lmod/linux-rhel7-ppc64le/Core"
      run "module load cuda python gcc/7.3.1 spectrum-mpi cmake/3.12.1 git"
      run "module load gsl numdiff random123 metis parmetis superlu-dist"
      run "module load trilinos netlib-lapack numdiff"
      run "module list"
      unset MPI_ROOT
      CXX=`which g++`
      CC=`which gcc`
      FC=`which gfortran`
    }
    ;;

  #------------------------------------------------------------------------------#
esac

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
