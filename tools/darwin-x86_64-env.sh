#!/bin/bash
#-----------------------------------------------------------------------------#
# Darwin Environment setups (Power9)
#-----------------------------------------------------------------------------#

source $draco_script_dir/darwin-env.sh

# symlinks will be generated for each machine that point to the correct
# installation directory.
export siblings="darwin-x86_64"

# The following toolchains will be used when releasing code:
environments="x86gcc730env x86intel1905env"

# Special setup for CTS-1: replace the 'draco-latest' symlink
(cd /usr/projects/$package; if [[ -L draco-latest ]]; then rm draco-latest; fi; ln -s $source_prefix draco-latest)

#--------------------------------------------------------------------------------------------------#
# Specify environments (modules)
#--------------------------------------------------------------------------------------------------#

case $ddir in

  #---------------------------------------------------------------------------#
  draco-7_7*|draco-7_8_*)
    function x86gcc730env()
    {
      export darwin_queue="-p volta-x86"
      run "module purge"
      module use --append ${VENDOR_DIR}/user_contrib
      module use --append ${VENDOR_DIR}-ec/Modules/$DRACO_ARCH
      module load user_contrib

      cflavor="gcc-7.3.0"
      mflavor="$cflavor-openmpi-3.1.3"
      lflavor="lapack-3.8.0"

      module_for_comp=`echo $cflavor | sed -e 's%-%/%'`
      module_for_mpi=`echo $mflavor | sed -e 's/.*open/open/' | sed -e 's%-%/%'`-`echo $cflavor | sed -e 's/-/_/'`
      noflavor="git ${module_for_comp} cmake/3.17.0"
      compflavor="gsl/2.5-$cflaavor netlib-lapack/3.8.0-$cflavor numdiff/5.9.0-$cflavor random123/1.09-$cflavor metis/5.1.0-$cflavor eospac/6.4.0-$cflavor ${module_for_mpi}"
      if [[ ${SLURM_JOB_PARTITION} =~ "volta-" ]] ; then
        # If we have GPUs on this node, then load the cuda toolkit module by default
        compflavor+=" cuda/10.2"
      fi
      mpiflavor="parmetis/4.0.3-$mflavor superlu-dist/5.2.2-$mflavor-${lflavor} trilinos/12.14.1-$mflavor-${lflavor} libquo/1.3-$mflavor"
      ec_mf="ndi csk/0.5.0-$cflavor"

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

    function x86intel1905env()
    {
      export darwin_queue="-p volta-x86"
      run "module purge"
      module use --append ${VENDOR_DIR}/user_contrib
      module use --append ${VENDOR_DIR}-ec/Modules/$DRACO_ARCH
      run "module load user_contrib gcc/7.3.0"
      unset CC CXX CPP FC F77 F90

      cflavor="intel-19.0.5"
      mflavor="$cflavor-openmpi-3.1.4"
      lflavor="lapack-3.8.0"

      module_for_comp=`echo $cflavor | sed -e 's%-%/%'`
      module_for_mpi=`echo $mflavor | sed -e 's/.*open/open/' | sed -e 's%-%/%'`-`echo $cflavor | sed -e 's/-/_/'`
      noflavor="git ${module_for_comp} cmake/3.17.0"
      compflavor="gsl/2.5-$cflavor netlib-lapack/3.8.0-$cflavor numdiff/5.9.0-$cflavor random123/1.09-$cflavor metis/5.1.0-$cflavor eospac/6.4.0-$cflavor ${module_for_mpi}"
      if [[ ${SLURM_JOB_PARTITION} =~ "volta-" ]] ; then
        # If we have GPUs on this node, then load the cuda toolkit module by default
        compflavor+=" cuda/10.2"
      fi
      mpiflavor="parmetis/4.0.3-$mflavor superlu-dist/5.2.2-$mflavor-${lflavor} trilinos/12.14.1-$mflavor-${lflavor} libquo/1.3-$mflavor"
      ec_mf="ndi csk/0.5.0-$cflavor"

      export dracomodules="$noflavor $compflavor $mpiflavor $ec_mf"
      for m in $dracomodules; do
        module load $m
      done
      # export CXX=`which g++`
      # export CC=`which gcc`
      # export FC=`which gfortran`
      export MPIEXEC_EXECUTABLE=`which mpirun`
      unset MPI_ROOT
      run "module list"
    }
    ;;

  *)
    die "darwin-x86_64-env.sh:: did not set any build environments, ddir = $ddir."
    ;;


  #---------------------------------------------------------------------------#

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

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
