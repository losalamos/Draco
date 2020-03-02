#!/bin/bash
#-----------------------------------------------------------------------------#
# Darwin Environment setups (Power9)
#-----------------------------------------------------------------------------#

source $draco_script_dir/darwin-env.sh

# symlinks will be generated for each machine that point to the correct
# installation directory.
export siblings="darwin-x86_64"

# The following toolchains will be used when releasing code:
environments="x86gcc730env"

# Special setup for CTS-1: replace the 'draco-latest' symlink
(cd /usr/projects/$package; if [[ -L draco-latest ]]; then rm draco-latest; fi; ln -s $source_prefix draco-latest)

#------------------------------------------------------------------------------#
# Specify environments (modules)
#------------------------------------------------------------------------------#

case $ddir in

  #---------------------------------------------------------------------------#
  draco-7_2* | draco-7_3* | draco-7_4*)
    function x86gcc730env()
    {
      export darwin_queue="-p volta-v100-x86"
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
    ;;

  *)
    die "darwin-x86_64-env.sh:: did not set any build environments, ddir = $ddir."
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
