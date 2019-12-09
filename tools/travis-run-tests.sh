#!/bin/bash -l
##---------------------------------------------------------------------------##
## File  : ./travis-run-tests.sh
## Date  : Tuesday, Jan 17, 2017, 15:55 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2017-2019, Triad National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# .travis.yml calls this script to build draco and run the tests.

# preliminaries and environment
set -e

cd ${SOURCE_DIR:-/home/travis/Draco}
source tools/common.sh

if [[ ${STYLE} ]]; then
  echo "checking style conformance..."

  # Return integer > 0 if 'develop' branch is found.
  function find_dev_branch
  {
    set -f
    git branch -a | grep -c develop
    set +f
  }

  # Ensure the 'develop' branch is available.  In some cases (merge a branch
  # that lives at github.com/lanl), the develop branch is missing in the
  # travis checkout. Since we only test files that are modified when comapred to
  # the 'develop' branch, the develop branch must be available locally.
  num_dev_branches_found=`find_dev_branch`
  if [[ $num_dev_branches_found == 0 ]]; then
    echo "no develop branches found."
    # Register the develop branch in draco/.git/config
    run "git config --local remote.origin.fetch +refs/heads/develop:refs/remotes/origin/develop"
    # Download the meta-data for the 'develop' branch
    run "git fetch"
    # Create a local tracking branch
    run "git branch -t develop origin/develop"
  fi

  # clang-format is installed at /usr/bin.
  export PATH=$PATH:/usr/bin
  # extract the TPL list from the Dockerfile
  export CLANG_FORMAT_VER="`grep \"ENV CLANG_FORMAT_VER\" tools/Dockerfile | sed -e 's/.*=//' | sed -e 's/\"//g'`"
  tools/check_style.sh -t

else
  echo "checking build and test..."

  # extract the TPL list from the Dockerfile
  export DRACO_TPL="`grep \"ENV DRACO_TPL\" tools/Dockerfile | sed -e 's/.*=//' | sed -e 's/\"//g'`"

  # Environment setup for the build...
  for item in ${DRACO_TPL}; do
    run "spack load ${item}"
  done

  if [[ $GCCVER ]]; then
    export CXX=`which g++-${GCCVER}`
    export CC=`which gcc-${GCCVER}`
    export FC=`which gfortran-${GCCVER}`
    export GCOV=`which gcov-${GCCVER}`
  else
    export CXX=`which g++`
    export CC=`which gcc`
    export FC=`which gfortran`
    export GCOV=`which gcov`
  fi
  echo "GCCVER = ${GCCVER}"
  echo "CXX    = ${CXX}"
  echo "FC     = ${FC}"
#  ls -1 /usr/bin/g*
#  if ! [[ $FC ]]; then export FC=gfortran; fi
  # ${FC} --version

  export OMP_NUM_THREADS=2
  if [[ ${WERROR} ]]; then
    for i in C CXX Fortran; do
      eval export ${i}_FLAGS+=\" -Werror\"
    done
  fi
  if [[ ${COVERAGE:-OFF} == "ON" ]]; then
    for i in C CXX Fortran; do
      eval export ${i}_FLAGS+=\" --coverage\"
    done
  fi

  # echo " "
  # echo "========== printenv =========="
  # printenv
  # echo " "

  if ! [[ $BUILD_DIR ]]; then die "BUILD_DIR not set by environment."; fi
  run "mkdir -p ${BUILD_DIR}"
  run "cd ${BUILD_DIR}"

  echo " "
  if [[ -f CMakeCache.txt ]]; then
    echo "===== CMakeCache.txt ====="
    run "cat CMakeCache.txt"
  fi
  echo "========"
  run "cmake -DDRACO_C4=${DRACO_C4} ${SOURCE_DIR}"
  echo "========"
  if [[ "${AUTODOC}" == "ON" ]]; then
    run "make autodoc"
    echo "========"
  else
    run "make -j 2"
    echo "========"
    # tstOMP_2 needs too many ppr (threads * cores) for Travis.
    run "ctest -j 2 -E \(c4_tstOMP_2\|c4_tstTermination_Detector_2\) --output-on-failure"
  fi
  cd -
  if [[ ${COVERAGE} == "ON" ]]; then
    echo "========"
    #which codecov
    #run "codecov --gcov-exec $GCOV"
    # https://docs.codecov.io/docs/testing-with-docker
    /bin/bash <(curl -s https://codecov.io/bash)
  fi
  echo "======== end .travis-run-tests.sh =========="
fi

# Finish up and report
