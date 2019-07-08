#!/bin/bash
#-----------------------------------------------------------------------------#
# ATS-2 Environment setups
#-----------------------------------------------------------------------------#

case $ddir in

  #---------------------------------------------------------------------------#
  draco-7_2_0)
    function p9gcc730env()
    {
      export darwin_queue="-p power9-asc -A asc-priority"
      export VENDOR_DIR=/projects/draco/vendors
      export DRACO_ARCH=`/usr/projects/draco/vendors/bin/target_arch`
      run "module purge"
      module use --append ${VENDOR_DIR}/user_contrib
      module use --append ${VENDOR_DIR}-ec/Modules/$DRACO_ARCH
      module load user_contrib


      cflavor="gcc-7.3.0"
      mflavor="$cflavor-openmpi-3.1.3"
      lflavor="lapack-3.8.0"
      noflavor="git gcc/7.3.0 cuda/10.1"
      compflavor="eospac/6.4.0-$cflavor cmake/3.14.2-$cflavor random123
numdiff gsl/2.5-$cflavor netlib-lapack/3.8.0-$cflavor metis/5.1.0-$cflavor
openmpi/p9/3.1.3-gcc_7.3.0"
      mpiflavor="parmetis/4.0.3-$mflavor
superlu-dist/5.2.2-${mflavor}-$lflavor trilinos/12.14.1-cuda-10.1-${mflavor}-$lflavor"

      # These aren't built for power architectures?
      # ec_mf="ndi eospac/6.3.0"

      # work around for known openmpi issues:
      # https://rtt.lanl.gov/redmine/issues/1229
      # eliminates warnings: "there are more than one active ports on host"
      # export OMPI_MCA_btl=^openib
      export UCX_NET_DEVICES=mlx5_0:1
      export UCX_WARN_UNUSED_ENV_VARS=n
      export OMPI_MCA_pml=ob1
      export OMPI_MCA_btl=self,vader

      export dracomodules="$noflavor $compflavor $mpiflavor $ec_mf"

      for m in $dracomodules; do
        module load $m
      done
      export CXX=`which g++`
      export CC=`which gcc`
      export FC=`which gfortran`
      export MPIEXEC_EXECUTABLE=`which mpirun`
      unset MPI_ROOT
      run "module list"
    }

    #-------------------------------------------------------------------------#

    function x86gcc730env()
    {
      export darwin_queue="-p volta-v100-x86"
      export VENDOR_DIR=/projects/draco/vendors
      export DRACO_ARCH=`/usr/projects/draco/vendors/bin/target_arch`
      run "module purge"
      module use --append ${VENDOR_DIR}/user_contrib
      module use --append ${VENDOR_DIR}-ec/Modules/$DRACO_ARCH
      module load user_contrib

      cflavor="gcc-7.3.0"
      mflavor="$cflavor-openmpi-3.1.3"
      lapackflavor="lapack-3.8.0"
      noflavor="emacs git ack gcc/7.3.0 cuda/10.1"
      compflavor="cmake/3.14.2-$cflavor gsl/2.5-$cflavor
netlib-lapack/3.8.0-$cflavor numdiff/5.9.0-$cflavor random123/1.09-$cflavor
metis/5.1.0-$cflavor eospac/6.4.0-$cflavor openmpi/3.1.3-gcc_7.3.0"
      mpiflavor="parmetis/4.0.3-$mflavor superlu-dist/5.2.2-$mflavor-$lapackflavor trilinos/12.14.1-cuda-10.1-$mflavor-$lapackflavor"
      ec_mf="ndi"

      export dracomodules="$noflavor $compflavor $mpiflavor $ec_mf"
      for m in $dracomodules; do
        module load $m
      done
      export CXX=`which g++`
      export CC=`which gcc`
      export FC=`which gfortran`
      export MPIEXEC_EXECUTABLE=`which mpirun`
      unset MPI_ROOT
      run "module list"
    }

    #-------------------------------------------------------------------------#

    function armgcc820env()
    {
      export darwin_queue="-p arm"
      export VENDOR_DIR=/projects/draco/vendors
      export DRACO_ARCH=`/usr/projects/draco/vendors/bin/target_arch`
      run "module purge"
      module use --append ${VENDOR_DIR}/user_contrib
      module use --append ${VENDOR_DIR}-ec/Modules/$DRACO_ARCH
      module load user_contrib

      cflavor="gcc-8.2.0"
      mflavor="$cflavor-openmpi-3.1.3"
      lapackflavor="lapack-3.8.0"
      noflavor="git gcc/8.2.0"
      compflavor="cmake/3.14.2-$cflavor gsl/2.5-$cflavor netlib-lapack/3.8.0-$cflavor numdiff/5.9.0-$cflavor random123/1.09-$cflavor metis/5.1.0-$cflavor eospac/6.4.0-$cflavor openmpi/3.1.3-gcc_8.2.0"
      mpiflavor="parmetis/4.0.3-$mflavor superlu-dist/5.2.2-$mflavor-$lapackflavor trilinos/12.14.1-$mflavor-$lapackflavor"
      ec_mf="ndi"

      export dracomodules="$noflavor $compflavor $mpiflavor $ec_mf"
      for m in $dracomodules; do
        module load $m
      done
      export CXX=`which g++`
      export CC=`which gcc`
      export FC=`which gfortran`
      export MPIEXEC_EXECUTABLE=`which mpirun`
      unset MPI_ROOT
      run "module list"
      # work around for known openmpi issues: https://rtt.lanl.gov/redmine/issues/1229
      export OMPI_MCA_btl=^openib
      export UCX_NET_DEVICES=mlx5_0:1

    }
    ;;

  #---------------------------------------------------------------------------#

esac

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
