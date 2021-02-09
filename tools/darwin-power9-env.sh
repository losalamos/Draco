#!/bin/bash
#--------------------------------------------------------------------------------------------------#
# Darwin Environment setups (Power9)
#--------------------------------------------------------------------------------------------------#

source $draco_script_dir/darwin-env.sh

# symlinks will be generated for each machine that point to the correct installation directory.
export siblings="darwin-power9"

# The following toolchains will be used when releasing code:
environments="p9gcc730env p9xl16117env"
#environments="p9gcc730env"

#--------------------------------------------------------------------------------------------------#
# Specify environments (modules)
#--------------------------------------------------------------------------------------------------#

case $ddir in

  #---------------------------------------------------------------------------#
  draco-7_9*)
    function p9gcc730env()
    {
      export darwin_queue="-p power9-asc -A asc-priority"
      run "module purge"
      run "module use --append /projects/draco/Modules"
      run "module load draco/power9-gcc730"
      run "module list"

      # work around for known openmpi issues:
      # https://rtt.lanl.gov/redmine/issues/1229
      # eliminates warnings: "there are more than one active ports on host"
      # export OMPI_MCA_btl=^openib
      export UCX_NET_DEVICES=mlx5_0:1
      export UCX_WARN_UNUSED_ENV_VARS=n
      export OMPI_MCA_pml=ob1
      export OMPI_MCA_btl=self,vader

      export CXX=`which g++`
      export CC=`which gcc`
      export FC=`which gfortran`
      export MPIEXEC_EXECUTABLE=`which mpirun`
      unset MPI_ROOT
    }

    function p9xl16117env()
    {
      export darwin_queue="-p power9-asc -A asc-priority"
      run "module purge"
      echo "VENDOR_DIR = $VENDOR_DIR"
      echo "DRACO_ARCH = $DRACO_ARCH"
      run "module use --append /projects/draco/Modules"
      run "module load draco/power9-xl16117"
      run "module list"

      # work around for known openmpi issues:
      # https://rtt.lanl.gov/redmine/issues/1229
      # eliminates warnings: "there are more than one active ports on host"
      # export OMPI_MCA_btl=^openib
      export UCX_NET_DEVICES=mlx5_0:1
      export UCX_WARN_UNUSED_ENV_VARS=n
      export OMPI_MCA_pml=ob1
      export OMPI_MCA_btl=self,vader
      export MPIEXEC_EXECUTABLE=`which mpirun`
      unset MPI_ROOT
    }
    ;;

  *) die "darwin-power9-env.sh:: did not set any build environments, ddir = $ddir." ;;
esac

#--------------------------------------------------------------------------------------------------#
# Sanity check
#--------------------------------------------------------------------------------------------------#

for env in $environments; do
  if [[ `fn_exists $env` -gt 0 ]]; then
    if [[ $verbose ]]; then echo "export -f $env"; fi
    export -f $env
  else
    die "Requested environment $env is not defined."
  fi
done

#--------------------------------------------------------------------------------------------------#
# End darwin-power9-env.sh
#--------------------------------------------------------------------------------------------------#
