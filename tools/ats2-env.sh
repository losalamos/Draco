#!/bin/bash
#------------------------------------------------------------------------------#
# ATS-2 Environment setups
#------------------------------------------------------------------------------#

export VENDOR_DIR=/usr/gapps/jayenne/vendors

# Sanity Check (Cray machines have very fragile module systems!)
if [[ -d $ParMETIS_ROOT_DIR ]]; then
  echo "ERROR: This script should be run from a clean environment."
  echo "       Try running 'rmdracoenv'."
  exit 1
fi

# symlinks will be generated for each machine that point to the correct
# installation directory.
if [[ `df | grep -c rz_gapps` -gt 0 ]]; then
  export siblings="rzansel"
else
  export siblings="sierra"
fi

# The following toolchains will be used when releasing code
environments="gcc731env xl20190820"

# Extra cmake options
export CONFIG_BASE+=" -DCMAKE_VERBOSE_MAKEFILE=ON"

# job launch options
case $siblings in
  rzansel) job_launch_options="-q pdebug" ;;
  sierra) job_launch_options="-q pbatch" ;;
  *) die "*FATAL ERROR* (ats2-env.sh) I only know how to set job_launch_options for rzansel and sierra." ;;
esac
export job_launch_options

# Special setup for ATS-1: replace the 'latest' symlink
(cd /usr/gapps/jayenne; if [[ -L draco-latest ]]; then rm draco-latest; fi; ln -s $source_prefix draco-latest)

#------------------------------------------------------------------------------#
# Specify environments (modules)
#------------------------------------------------------------------------------#

if ! [[ $ddir ]] ;then
  echo "FATAL ERROR: Expected ddir to be set in the environment. (ats2-env.sh)"
  exit 1
fi

case $ddir in

  #------------------------------------------------------------------------------#
  draco-7_4*)
    function gcc731env()
    {
      run "module purge"
      run "module use /usr/gapps/jayenne/vendors-ec/spack.20190616/share/spack/lmod/linux-rhel7-ppc64le/Core"
      run "module load StdEnv"
      run "module unload spectrum-mpi xl"
      run "module load gcc/7.3.1 spectrum-mpi/2019.06.24"
      run "module load python/3.7.2 cmake/3.14.5 git gsl numdiff random123 metis netlib-lapack"
      run "module load parmetis superlu-dist trilinos csk"
      run "module load eospac/6.4.0"
      # ndi quo
      run "module list"
      unset MPI_ROOT
      CXX=`which g++`
      CC=`which gcc`
      FC=`which gfortran`
    }
    function xl20190820()
    {
      run "module purge"
      run "module use /usr/gapps/jayenne/vendors-ec/spack.20190616/share/spack/lmod/linux-rhel7-ppc64le/Core"
      run "module load StdEnv"
      run "module unload spectrum-mpi"
      run "module load spectrum-mpi/2019.06.24"
      run "module load python/3.7.2 cmake/3.14.5 git gsl numdiff random123 metis netlib-lapack"
      run "module load parmetis superlu-dist csk"
      run "module load eospac/6.4.0"
      # ndi quo
      run "module list"
      unset MPI_ROOT
    }
    ;;

  #------------------------------------------------------------------------------#
  draco-7_2* | draco-7_3* )
    function gcc731env()
    {
      run "module purge"
      run "module use /usr/gapps/jayenne/vendors-ec/spack.20190616/share/spack/lmod/linux-rhel7-ppc64le/Core"
      run "module load StdEnv cuda python/3.7.2 gcc/7.3.1 spectrum-mpi/2019.01.30"
      run "module load cmake/3.14.5 git gsl numdiff random123 metis netlib-lapack"
      run "module load parmetis superlu-dist trilinos csk"
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
  draco-7_1*)
    function gcc731env()
    {
      run "module purge"
      run "module use /usr/gapps/user_contrib/spack.20190314/share/spack/lmod/linux-rhel7-ppc64le/Core"
      run "module load StdEnv cuda python gcc/7.3.1 spectrum-mpi cmake/3.12.1 git"
      run "module load gsl numdiff random123 metis parmetis superlu-dist"
      run "module load trilinos netlib-lapack numdiff"
      run "module list"
      unset MPI_ROOT
      CXX=`which g++`
      CC=`which gcc`
      FC=`which gfortran`
    }
    ;;

  *)
    die "ats2-env.sh:: did not set any build environments, ddir = $ddir."
    ;;

  #------------------------------------------------------------------------------#
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
