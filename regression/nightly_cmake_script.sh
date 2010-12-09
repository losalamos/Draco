#!/bin/bash

# The work_dir is the location for the source and build directories
# /home/regress/draco/cmake_draco/
#      source/  <-- Source files checked out from CVS go here.
#      build/   <-- Make is run in this location.

unset http_proxy

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

export work_dir=/home/regress/cmake_draco/Nightly/Release
ctest -VV -S ~/draco/regression/Draco_gcc.cmake,${dashboard_type},Release

# # Debug build
export work_dir=/home/regress/cmake_draco/Nightly/Debug
ctest -VV -S ~/draco/regression/Draco_gcc.cmake,${dashboard_type},Debug

# Coverage build
module load bullseyecoverage
export work_dir=/home/regress/cmake_draco/Nightly/Coverage
ctest -VV -S ~/draco/regression/Draco_gcc.cmake,${dashboard_type},Debug,Coverage
module unload bullseyecoverage