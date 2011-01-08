#!/bin/bash

# The work_dir is the location for the source and build directories
# /home/regress/draco/cmake_draco/
#      source/  <-- Source files checked out from CVS go here.
#      build/   <-- Make is run in this location.

# Environment setup

unset http_proxy
export VENDOR_DIR=/ccs/codes/radtran/vendors/Linux64

if test -z "$MODULESHOME"; then
  # This is a new login
  if test -f /home/regress/environment/Modules/init/bash; then
    source /home/regress/environment/Modules/init/bash
    module load grace BLACS SCALAPACK SuperLU_DIST/2.4 gandolf gcc/4.3.4 gsl hypre
    module load lapack ndi openmpi ParMetis/3.1.1 trilinos/10.4.0 cmake
    module list
  fi
fi


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

dashboard_type=Nightly
base_dir=/home/regress/cmake_draco
script_dir=/home/regress/cmake_draco

# Release build
build_type=Release
export work_dir=${base_dir}/${dashboard_type}/${build_type}
ctest -VV -S ${script_dir}/regression/Draco_gcc.cmake,${dashboard_type},${build_type}
(cd $work_dir/build; make install)


# Debug build
build_type=Debug
export work_dir=${base_dir}/${dashboard_type}/${build_type}
ctest -VV -S ${script_dir}/regression/Draco_gcc.cmake,${dashboard_type},${build_type}
(cd $work_dir/build; make install)

# Coverage build
module load bullseyecoverage
CXX=`which g++`
CC=`which gcc`
export work_dir=${base_dir}/${dashboard_type}/Coverage
ctest -VV -S ${script_dir}/regression/Draco_gcc.cmake,${dashboard_type},Debug,Coverage
module unload bullseyecoverage
(cd $work_dir/build; make install)
