#!/bin/bash
#------------------------------------------------------------------------------#
# CTS-1 Environment setups
#------------------------------------------------------------------------------#

export VENDOR_DIR=/usr/projects/draco/vendors

# symlinks will be generated for each machine that point to the correct
# installation directory.
if [[ $(df | grep yellow | grep -c jayenne) -gt 0 ]]; then
  export siblings="snow badger kodiak grizzly"
else
  export siblings="fire ice cyclone"
fi

# The following toolchains will be used when releasing code
environments="intel1904env intel1704env gcc830env gcc740env"

# Extra cmake options
export CONFIG_BASE+=" -DCMAKE_VERBOSE_MAKEFILE=ON"

# SLURM
avail_queues=$(sacctmgr -np list assoc user="$LOGNAME" | sed -e 's/.*|\(.*dev.*\|.*access.*\)|.*/\1/' | sed -e 's/|.*//')
case "$avail_queues" in
  *access*) access_queue="-A access --qos=access" ;;
  *dev*) access_queue="--qos=dev" ;;
esac
export access_queue

# Special setup for CTS-1: replace the 'latest' symlink
if [[ ${package:-false} == false ]] ; then die "package not defined"; fi
if [[ ${source_prefix:-false} == false ]] ; then die "source_prefix not defined"; fi
(cd "/usr/projects/${package:=draco}" || exit; if [[ -L latest ]]; then rm latest; fi; ln -s "${source_prefix:-source_prefix}" latest)

#------------------------------------------------------------------------------#
# Specify environments (modules)
#------------------------------------------------------------------------------#

if [[ ${ddir:=false} == false ]] ;then
  echo "FATAL ERROR: Expected ddir to be set in the environment. (cts1-env.sh)"
  exit 1
fi

case "${ddir}" in

  #------------------------------------------------------------------------------#
  draco-7_5*|draco-7_6*|draco-7_7*)
    function intel1904env
    {
      run "module purge"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load friendly-testing user_contrib"
      run "module load cmake/3.17.0 git numdiff python/3.6-anaconda-5.0.1"
      run "module load intel/19.0.4 openmpi/2.1.2"
      run "unset MPI_ROOT"
      run "module load random123 eospac/6.4.0 gsl"
      run "module load mkl metis ndi csk qt quo"
      run "module load parmetis superlu-dist/5.1.3 trilinos/12.10.1"
      run "module list"
    }
    function intel1704env
    {
      run "module purge"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load friendly-testing user_contrib"
      run "module load cmake/3.17.0 git numdiff python/3.6-anaconda-5.0.1"
      run "module load intel/17.0.4 openmpi/2.1.2"
      run "unset MPI_ROOT"
      run "module load random123 eospac/6.4.0 gsl"
      run "module load mkl metis ndi csk qt quo"
      run "module load parmetis superlu-dist/5.1.3 trilinos/12.10.1"
      run "module list"
    }
    function gcc830env()
    {
      run "module purge"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load friendly-testing user_contrib"
      run "module load cmake/3.17.0 git numdiff python/3.6-anaconda-5.0.1"
      run "module load gcc/8.3.0 openmpi/2.1.2"
      run "unset MPI_ROOT"
      run "module load random123 eospac/6.4.0 gsl"
      run "module load mkl metis ndi csk qt quo"
      run "module load parmetis superlu-dist/5.1.3 trilinos/12.10.1"
      run "module list"
    }
    function gcc740env()
    {
      run "module purge"
      run "module use --append ${VENDOR_DIR}-ec/modulefiles"
      run "module load friendly-testing user_contrib"
      run "module load cmake/3.17.2 git numdiff python/3.6-anaconda-5.0.1"
      run "module load gcc/7.4.0 openmpi/2.1.2"
      run "unset MPI_ROOT"
      run "module load random123 eospac/6.4.0 gsl"
      run "module load mkl metis ndi csk qt quo"
      run "module load parmetis superlu-dist/5.1.3 trilinos/12.10.1"
      run "module list"
    }
    ;;

#------------------------------------------------------------------------------#

  *)
    die "cts1-env.sh:: did not set any build environments, ddir = $ddir."
    ;;

esac

#------------------------------------------------------------------------------#
# Sanity check
#------------------------------------------------------------------------------#

for env in ${environments}; do
  if [[ $(fn_exists "$env") -gt 0 ]]; then
    if [[ "${verbose:-false}" != false ]]; then echo "export -f $env"; fi

    export -f "${env?}"
  else
    die "Requested environment $env is not defined."
  fi
done

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
