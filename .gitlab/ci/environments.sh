#!/bin/bash -l
# -*- Mode: sh -*-
#--------------------------------------------------------------------------------------------------#
# File  : .gitlab/ci/environments.sh
# Date  : Monday, Jun 01, 2020, 15:43 pm
# Author: Kelly Thompson
# Note  : Copyright (C) 2020-2021, Triad National Security, LLC., All rights are reserved.
#--------------------------------------------------------------------------------------------------#

echo "==> Setting up CI environment..."
echo "    SITE_ID = ${SITE_ID}"
[[ -n "${SCHEDULER_PARAMETERS}" ]] && echo "    Using: salloc ${SCHEDULER_PARAMETERS}"

case ${SITE_ID} in
  darwin | ccscs*) ;;
  *) die ".gitlab/ci/environments.sh :: SITE_ID not recognized, SITE_ID = ${SITE_ID}" ;;
esac

echo "    DRACO_ENV = ${DRACO_ENV}"
echo "    ARCH     = ${ARCH}"

if [[ "${SITE_ID}" == "darwin" ]]; then
  DRACO_ARCH=$(/usr/projects/draco/vendors/bin/target_arch)
  run "module use --append /projects/draco/Modules"
  case ${DRACO_ENV} in
    arm-gcc930 | power9-gcc930-smpi | power9-xl16117 | x64-gcc930 | x64-intel1905 ) ;;
    *) die ".gitlab/ci/environments.sh :: DRACO_ENV not recognized, DRACO_ENV = ${DRACO_ENV}" ;;
  esac
  if [[ "${MPIARCH:-notset}" == "openmpi" ]]; then
    disable_openib=$(sinfo -N -n "${HOSTNAME}" -o %all | grep -c ib:none)
    if [[ ${disable_openib} != 0 ]]; then
      export MPI_PREFLAGS="--mca btl ^openib"
    fi
  fi
  export DRACO_ARCH

elif [[ "${SITE_ID}" =~ "ccscs" ]]; then
  run "module use --append /ccs/codes/radtran/Modules"
  export PATH=/scratch/vendors/bin:$PATH # clang-format
  case ${DRACO_ENV} in
    gcc1020* | llvm11 ) ;;
    *) die ".gitlab/ci/environments.sh :: DRACO_ENV not recognized, DRACO_ENV = ${DRACO_ENV}" ;;
  esac
fi

#--------------------------------------------------------------------------------------------------#
# Setup requested CI environment
#--------------------------------------------------------------------------------------------------#

run "module load draco/${DRACO_ENV}"
if [[ "${SLURM_JOB_PARTITION}" =~ "volta" || "${SLURM_JOB_PARTITION}" =~ "gpu" ]]; then
  module load cuda/11.0
fi
run "module list"

#--------------------------------------------------------------------------------------------------#
# End environments.sh
#--------------------------------------------------------------------------------------------------#
