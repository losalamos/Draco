#!/bin/bash
#--------------------------------------------------------------------------------------------------#
# Darwin Environment setups (x86_64 + gpu)
#--------------------------------------------------------------------------------------------------#

source $draco_script_dir/darwin-env.sh

# symlinks will be generated for each machine that point to the correct installation directory.
export siblings="darwin-x86_64"

# The following toolchains will be used when releasing code:
environments="x86gcc930env x86intel1905env"

# Special setup for CTS-1: replace the 'draco-latest' symlink
(cd /usr/projects/$package; if [[ -L draco-latest ]]; then rm draco-latest; fi; \
ln -s $source_prefix draco-latest)

#--------------------------------------------------------------------------------------------------#
# Specify environments (modules)
#--------------------------------------------------------------------------------------------------#

case $ddir in

  #---------------------------------------------------------------------------#
  draco-7_9*)
    function x86gcc930env()
    {
      export darwin_queue="-p volta-x86"
      run "module purge"
      run "module use --append /projects/draco/Modules"
      run "module load draco/x64-gcc930"
      if [[ ${SLURM_JOB_PARTITION} =~ "volta-" || ${SLURM_JOB_PARTITION} =~ "gpu" ]]; then
        run "module load cuda/11.2-beta"
      fi
      run "module list"
      export CXX=`which g++`
      export CC=`which gcc`
      export FC=`which gfortran`
      export MPIEXEC_EXECUTABLE=`which mpirun`
      unset MPI_ROOT
    }

    function x86intel1905env()
    {
      export darwin_queue="-p volta-x86"
      run "module purge"
      run "module use --append /projects/draco/Modules"
      run "module load draco/x64-intel1905"
      #if [[ ${SLURM_JOB_PARTITION} =~ "volta-" || ${SLURM_JOB_PARTITION} =~ "gpu" ]]; then
      #  run "module load cuda/11.2-beta"
      #fi
      run "module list"
      export MPIEXEC_EXECUTABLE=`which mpirun`
      unset MPI_ROOT
    }
    ;;

  *) die "darwin-x86_64-env.sh:: did not set any build environments, ddir = $ddir." ;;
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
# End darwin-x86_64-env.sh
#--------------------------------------------------------------------------------------------------#
