#!/bin/bash

# IMPORTANT:
# The cron job that starts this file should start no earlier than 1
# am.  This time is coordinated with the CTEST_NIGHTLY_START_TIME
# which is fixed wrt DST.  If this guidance is not followed, you may
# see svn updates tagged from 2 days ago.

# The work_dir is the location for the source and build directories
# /home/regress/draco/cmake_draco/
#      source/  <-- Source files checked out from SVN go here.
#      build/   <-- Make is run in this location.

# Environment setup

umask 0002
unset http_proxy
unset HTTP_PROXY
export VENDOR_DIR=/ccs/codes/radtran/vendors/Linux64

(cd /home/regress/environment/Modules; svn update; svn status -u)
(cd /home/regress/cmake_draco/regression; svn update; svn status -u)

if test -z "$MODULESHOME"; then
  # This is a new login
  if test -f /ccs/codes/radtran/vendors/modules-3.2.9/init/bash; then
    source /ccs/codes/radtran/vendors/modules-3.2.9/init/bash
  fi
fi
module load grace gcc gsl lapack openmpi cmake numdiff svn valgrind
module list

echo " "
echo "top -b -n 1"
top -b -n 1 | head -n 15

# Run the ctest (regression) script.  This script will take the following build steps: 
# 1. svn update
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
    base_dir=/home/regress/cmake_draco
    script_dir=/home/regress/cmake_draco
else
    dashboard_type=Experimental
    base_dir=/var/tmp/${USER}/regress/cmake_draco
    script_dir=/ccs/codes/radtran/regress-homedir/cmake_draco
fi

# compiler
comp=gcc

# Release build
build_type=Release
export work_dir=${base_dir}/${dashboard_type}_${comp}/${build_type}
mkdir -p ${work_dir}
ctest -VV -S ${script_dir}/regression/Draco_gcc.cmake,${dashboard_type},${build_type}

# Debug build
build_type=Debug
export work_dir=${base_dir}/${dashboard_type}_${comp}/${build_type}
mkdir -p ${work_dir}
ctest -VV -S ${script_dir}/regression/Draco_gcc.cmake,${dashboard_type},${build_type}

# Coverage build
build_type=Coverage
module load bullseyecoverage/8.7.17
CXX=`which g++`
CC=`which gcc`
export work_dir=${base_dir}/${dashboard_type}_${comp}/${build_type}
mkdir -p ${work_dir}
ctest -VV -S ${script_dir}/regression/Draco_gcc.cmake,${dashboard_type},Debug,${build_type}
module unload bullseyecoverage

