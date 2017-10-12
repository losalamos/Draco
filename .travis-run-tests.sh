#!/bin/bash
##---------------------------------------------------------------------------##
## File  : ./travis-run-tests.sh
## Date  : Tuesday, Jan 17, 2017, 15:55 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# .travis.yml calls this script to build draco and run the tests.

# preliminaries and environment
set -e
source regression/scripts/common.sh

topdir=`pwd` # /home/travis/build/lanl/Draco
# HOME = /home/travis
# USER = travis
# GROUP = travis
RANDOM123_VER=1.09
CMAKE_VERSION=3.9.0-Linux-x86_64
NUMDIFF_VER=5.8.1
CLANG_FORMAT_VER=3.9
OPENMPI_VER=1.10.5
GCCVER=6
export CXX=`which g++-${GCCVER}`
export CC=`which gcc-${GCCVER}`
export FC=`which gfortran-${GCCVER}`

if [[ ${STYLE} ]]; then
  regression/check_style.sh -t
else
  mkdir -p build
  cd build
  # configure
  # -Wl,--no-as-needed is a workaround for bugs.debian.org/457284 .
  echo " "
  echo "${CMAKE} -DCMAKE_EXE_LINKER_FLAGS=\"-Wl,--no-as-needed\" .."
  ${CMAKE} -DCMAKE_EXE_LINKER_FLAGS="-Wl,--no-as-needed" ..
  error_code=$?
  # if configure was successful, the start the build, otherwise abort.
  if [[ $error_code -eq 0 ]]; then
    echo " "
    echo "make -j 2 VERBOSE=1"
    make -j 2 VERBOSE=1
    error_code=$?
  else
    echo "configure failed, errorcode=$error_code"
    exit $error_code
  fi
  # if the build was successful, then run the tests, otherwise abort.
  if [[ $error_code -eq 0 ]]; then
    echo " "
    echo "${CTEST} -VV -E \(c4_tstOMP_2\)"
    ${CTEST} -VV -E \(c4_tstOMP_2\)
    error_code=$?
  else
    echo "build failed, errorcode=$error_code"
    exit $error_code
  fi
  # if some of the tests failed, rerun them with full verbosity to provide more
  # information in the output log.
  # if ! [[ $error_code -eq 0 ]]; then
  #   echo "some tests failed, errorcode=$error_code, trying again with more verbosity."
  #   echo " "
  #   echo "${CTEST} --rerun-failed -VV -E \(c4_tstOMP_2\)"
  #   ${CTEST} --rerun-failed -VV -E \(c4_tstOMP_2\)
  #   exit $error_code
  #   echo "error_code=$error_code"
  # fi
  echo " "
  echo "make install DESTDIR=${HOME}"
  make install DESTDIR="${HOME}"
  cd -
fi

# Finish up and report
