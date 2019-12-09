#!/bin/bash
#-----------------------------------------------------------------------------#
# Darwin Environment setups (ARM)
#-----------------------------------------------------------------------------#

source $draco_script_dir/darwin-env.sh

# symlinks will be generated for each machine that point to the correct
# installation directory.
export siblings="darwin-arm"

# The following toolchains will be used when releasing code:
environments="armgcc820env"

#------------------------------------------------------------------------------#
# Specify environments (modules)
#------------------------------------------------------------------------------#

case $ddir in

  #---------------------------------------------------------------------------#
  draco-7_2* | draco-7_3*)
    function armgcc820env()
    {
      export darwin_queue="-p arm"
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
      # ec_mf="ndi"

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

  *)
    die "darwin-arm-env.sh:: did not set any build environments, ddir = $ddir."
    ;;

  #---------------------------------------------------------------------------#

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
