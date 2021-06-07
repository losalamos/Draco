#!/bin/bash -l
#--------------------------------------------------------------------------------------------------#
# File  : ./.gitlab/ci/gitlab-ci-run-tests.sh
# Date  : Tuesday, Jun 02, 2020, 12:28 pm
# Author: Kelly Thompson <kgt@lanl.gov>
# Note  : Copyright (C) 2020-2021, Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#

# preliminaries and environment
set -e
# shellcheck source=.gitlab/ci/common.sh
source "${JAYENNE_SOURCE_DIR}/.gitlab/ci/common.sh"
# shellcheck source=.gitlab/ci/environments.sh
source "${JAYENNE_SOURCE_DIR}/.gitlab/ci/environments.sh"

if [[ $(which lscpu | grep -c lscpu) -gt 0 ]]; then
  # run "lscpu"
  # NPROC=`lscpu | grep CPU\(s\) | head -n 1 | awk '{ print $2 }'`
  NPROC=$(grep -c processor < /proc/cpuinfo)
fi

[[ -z "${CTEST_NPROC}" ]] && CTEST_NPROC=$NPROC || echo "limiting CTEST_NPROC = $CTEST_NPROC"
[[ -z "${MAXLOAD}" ]] && MAXLOAD=$NPROC || echo "limiting MAXLOAD = $MAXLOAD"

#if [[ -d /ccs/opt/texlive/2018/texmf-dist/tex/latex/newunicodechar ]]; then
#  TEXINPUTS="${TEXINPUTS:+${TEXINPUTS}:}/ccs/opt/texlive/2018/texmf-dist/tex/latex/newunicodechar"
#fi

# Report number of cores/threads, etc.
# run "lstopo --no-io --no-bridges --no-icaches"

echo -e "\n========== printenv ==========\n"
[[ -z "${SLURM_NODELIST}" ]] || echo "SLURM_NODELIST = ${SLURM_NODELIST}"
echo "HOSTNAME       = ${HOSTNAME}"
echo -e "NPROC       = ${NPROC}\n"
echo -e "CTEST_NPROC = ${CTEST_NPROC}\n"
echo -e "MAXLOAD     = ${MAXLOAD}\n"
echo -e "EXTRA_CMAKE_ARGS = ${EXTRA_CMAKE_ARGS}\n"
echo -e "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}\n"
echo -e "EXTRA_CTEST_ARGS = ${EXTRA_CTEST_ARGS}\n"
if [[ "${AUTODOC}" == "ON" ]]; then
  run "which doxygen"
  run "which latex"
  run "which pdflatex"
  echo "TEXINPUTS = ${TEXINPUTS}"
fi
if [[ "${EXTRA_CTEST_ARGS}" =~ memcheck ]]; then
  run "which valgrind"
fi
printenv >& environment.log
run "pwd"

[[ "${EXTRA_CMAKE_ARGS}" =~ "CODE_COVERAGE" ]] && CODECOV=ON || CODECOV=OFF

#--------------------------------------------------------------------------------------------------#
# Setup compiler flags
#--------------------------------------------------------------------------------------------------#

for i in C CXX Fortran CUDA; do

  # Turn build warnings into fatal errors
  if [[ "${DRACO_ENV}" =~ "xl16" ]]; then
    eval export ${i}_FLAGS+=\" -qhalt=w\"
  else
    if [[ "${i}" == "Fortran" ]] && [[ "${FC}" =~ "ifort" ]]; then
      # do not add the flag
      echo "Skip adding -Werror to ${i}_FLAGS (FC = ${FC})"
    elif [[ "${i}" == "CUDA" ]]; then
      eval export ${i}_FLAGS+=\" -Werror all-warnings\"
    else
      # Works for gcc, llvm, and icpc (not ifort).
      eval export ${i}_FLAGS+=\" -Werror\"
    fi
  fi

  # Enable coverage for Debug builds
  if [[ ${CMAKE_BUILD_TYPE} == Debug ]] && [[ "${CODECOV}" == "ON" ]] ; then
    eval export ${i}_FLAGS+=\" --coverage\"
  fi
done
echo "C_FLAGS       = $C_FLAGS"
echo "CXX_FLAGS     = $CXX_FLAGS"
# shellcheck disable=2154
echo "Fortran_FLAGS = $Fortran_FLAGS"
echo "CUDA_FLAGS = $CUDA_FLAGS"

# Tweak MPI settings
# export OMPI_MCA_btl=self,sm
export OMPI_MCA_btl=^openib

#--------------------------------------------------------------------------------------------------#
# Draco (build and install libraries only, no tests)
#--------------------------------------------------------------------------------------------------#

echo -e "\n========== Draco ==========\n"
if ! [[ -d "${DRACO_BINARY_DIR}" ]]; then run "mkdir -p ${DRACO_BINARY_DIR}"; fi
if ! [[ -d "${DRACO_SOURCE_DIR}" ]]; then
  run "cd ${CI_PROJECT_DIR} && git clone https://github.com/lanl/Draco draco"
else
  run "cd ${DRACO_SOURCE_DIR} && git pull"
fi
run "cd ${DRACO_BINARY_DIR}"
echo -e "\n========"
run "cmake ${EXTRA_CMAKE_ARGS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DBUILD_TESTING=OFF ${DRACO_SOURCE_DIR}"
echo " "
run "make -j -l $MAXLOAD"
# don't echo all of the install commands to the terminal.
# shellcheck disable=2086
echo "make -j -l $MAXLOAD install"
make -j -l "$MAXLOAD" install &> /dev/null

#--------------------------------------------------------------------------------------------------#
# Jayenne
#--------------------------------------------------------------------------------------------------#

echo -e "\n========== Jayenne ==========\n"
if ! [[ -d ${JAYENNE_BINARY_DIR} ]]; then run "mkdir -p ${JAYENNE_BINARY_DIR}"; fi
run "cd ${JAYENNE_BINARY_DIR}"
echo -e "\n========"
run "cmake ${EXTRA_CMAKE_ARGS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} ${JAYENNE_SOURCE_DIR}"
echo " "
#if [[ -f CMakeCache.txt ]]; then
#  echo -e "\n===== CMakeCache.txt =====\n"
#  run "cat CMakeCache.txt"
#fi
echo -e "\n========\n"
if [[ "${AUTODOC}" == "ON" ]]; then
  run "make -j -l $MAXLOAD autodoc"
else
  run "make -j -l $MAXLOAD" #  VERBOSE=1
  echo -e "\n========\n"
  run "ctest ${EXTRA_CTEST_ARGS} -j ${CTEST_NPROC} --test-load ${MAXLOAD} --output-on-failure --stop-on-failure"
fi

# Generate a coverage report (text and html)
if [[ "${CODECOV}" == "ON" ]]; then
  run "make covrep"
fi

if [[ ${CTEST_MODE} == "Nightly" ]]; then
  echo -e "\n========\n"
  echo "Nighlty regressions should run a ctest in script mode and post to CDash"
fi

echo -e "\n======== end .gitlab-ci-run-tests.sh ==========\n"

#--------------------------------------------------------------------------------------------------#
# End .gitlab/ci/gitlab-ci-run-tests.sh
#--------------------------------------------------------------------------------------------------#
