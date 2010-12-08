#!/bin/bash

# The work_dir is the location for the source and build directories
# /home/regress/draco/cmake_draco/
#      source/  <-- Source files checked out from CVS go here.
#      build/   <-- Make is run in this location.

unset http_proxy
export work_dir=/home/regress/draco/cmake_draco

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
ctest -VV -S $work_dir/source/regression/Draco_gcc.cmake,Nightly,Release

