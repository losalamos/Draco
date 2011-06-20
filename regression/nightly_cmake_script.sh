#!/bin/bash

# IMPORTANT:
# The cron job that starts this file should start no earlier than 1
# am.  This time is coordinated with the CTEST_NIGHTLY_START_TIME
# which is fixed wrt DST.  If this guidance is not followed, you may
# see cvs updates tagged from 2 days ago.

# The work_dir is the location for the source and build directories
# /home/regress/draco/cmake_draco/
#      source/  <-- Source files checked out from CVS go here.
#      build/   <-- Make is run in this location.

# Environment setup

umask 0002
unset http_proxy
export VENDOR_DIR=/ccs/codes/radtran/vendors/Linux64

(cd /home/regress/environment/Modules; cvs -q update -AdP)

if test -z "$MODULESHOME"; then
  # This is a new login
  if test -f /home/regress/environment/Modules/init/bash; then
    source /home/regress/environment/Modules/init/bash
    module load grace BLACS SCALAPACK SuperLU_DIST/2.4 gandolf gcc/4.3.4 gsl 
    module load lapack-gcc ndi openmpi ParMetis/3.1.1 trilinos/10.4.0 cmake
    module load valgrind hypre
    module list
  fi
fi
module load valgrind

# Run the ctest (regression) script.  This script will take the following build steps: 
# 1. cvs update
# 2. run cmake to build Makefiles
# 3. run make to build libraries and tests
# 4. Run the unit tests
# 5. Post the results to coder.lanl.gov/cdash
#
# Options are:
# Regression type: Experimental (default), Nightly, Continuous
# Build type     : Release, Debug

if test "`whoami`" == "regress"; then
    dashboard_type=Nightly
else
    dashboard_type=Experimental
fi

base_dir=/home/regress/cmake_draco
script_dir=/home/regress/cmake_draco

# compiler
comp=gcc

# Release build
build_type=Release
export work_dir=${base_dir}/${dashboard_type}_${comp}/${build_type}
ctest -VV -S ${script_dir}/regression/Draco_gcc.cmake,${dashboard_type},${build_type}

# Debug build
build_type=Debug
export work_dir=${base_dir}/${dashboard_type}_${comp}/${build_type}
ctest -VV -S ${script_dir}/regression/Draco_gcc.cmake,${dashboard_type},${build_type}

# Coverage build
build_type=Coverage
module load bullseyecoverage
CXX=`which g++`
CC=`which gcc`
export work_dir=${base_dir}/${dashboard_type}_${comp}/${build_type}
ctest -VV -S ${script_dir}/regression/Draco_gcc.cmake,${dashboard_type},Debug,${build_type}
module unload bullseyecoverage

