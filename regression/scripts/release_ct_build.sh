#!/bin/bash -l

#MSUB -l walltime=01:00:00
#MSUB -l nodes=1:ppn=16
#MSUB -j oe
#MSUB -o /usr/projects/draco/draco-6_13_0/logs/release_ct_build_intel14.log

#----------------------------------------------------------------------#
# The script starts here
#----------------------------------------------------------------------#

echo "Here we go..." > /usr/projects/draco/draco-6_13_0/logs/release_ct_build_intel14.log

# Permissions - new files should be marked u+rwx,g+rwx,o+rx
umask 0002
build_permissions="g+rwX"
install_permissions="g+rwX,o=g-w"

# environment (use draco modules)
dirs="/usr/projects/jayenne/regress/draco/environment/Modules/hpc \
/usr/projects/jayenne/regress/draco/environment/Modules/ct-fe"
for dir in $dirs; do 
  if test -z `echo $MODULEPATH | grep $dir`; then 
    echo "module use $dir"
    module use $dir
  fi
done

module unload PrgEnv-intel PrgEnv-pgi
module unload cmake numdiff svn gsl
module unload papi perftools
module load PrgEnv-intel
# module swap intel intel/13.1.3.192
module unload xt-libsci xt-totalview 
module load gsl/1.14 lapack/3.4.1
module load cmake numdiff subversion emacs
module load trilinos SuperLU_DIST
module load ParMetis ndi random123 eospac/v6.2.4beta.3
module list

export OMP_NUM_THREADS=8

# Define your source and build information here.

ddir="draco-6_13_0"
platform="ct"
dmpi=craympich2
df90=intel1402
dcpp=intel1402

source_prefix="/usr/projects/draco/$ddir"
build_prefix="/lscratch1/$USER/$ddir/$platform-${dmpi}-${df90}-${dcpp}"
install_prefix="$source_prefix/$platform-${dmpi}-${df90}-${dcpp}"

MAKE="make -j48"
CONFIG_BASE="-DDRACO_VERSION_PATCH=0"
CC=`which cc`
CXX=`which CC`
FC=`which ftn`
printenv

# =============================================================================
# checkit_core
#
# This script can be used either stand-alone or sourced into a buildit
# script which defines the variables that it needs.
#
# It requires the following enviromnent variables to be defined:
#
# |----------------+-------------+--------------------------------------------|
# | Name           | Specificity | Description                                |
# |----------------+-------------+--------------------------------------------|
# | build_prefix   | build       | Destination directory of the build trees.  |
# | MAKE           | platform    | Make command.                              |
# |----------------+-------------+--------------------------------------------|
# 
# Optionally, if "dry_run" is defined, all commands are output to
# stdio, but not executed. This can be useful to create a "rebuild"
# script which repeats the configure and build process. 
#
#
# Results:
#
# This script runs $MAKE run in four versions of the code: debug,
# debug_nodbc, opt, and opt_log.
#
# =============================================================================

# Helpful functions:
die () { echo "ERROR: $1"; exit 1;}

run () {
    echo $1
    if ! test $dry_run; then
    eval $1
    fi
}

# Define the meanings of various configure features:
# DBC_OFF="-DDRACO_DBC_LEVEL=0"
# DBC_ON="-DDRACO_DBC_LEVEL=7"

OPTIMIZE_ON="-DCMAKE_BUILD_TYPE=Release"
OPTIMIZE_OFF="-DCMAKE_BUILD_TYPE=Debug"

LOGGING_ON="-DDRACO_DIAGNOSTICS=7 -DDRACO_TIMING=1"
LOGGING_OFF="-DDRACO_DIAGNOSTICS=0 -DDRACO_TIMING=0"

# Define the meanings of the various code versions:

VERSIONS=( "debug" "opt" )
OPTIONS=(\
    "$OPTIMIZE_OFF $LOGGING_OFF" \
    "$OPTIMIZE_ON  $LOGGING_OFF" \
)

PACKAGES=("draco")

# =============
# Configure, Build and Run the Tests
# =============


# Loop over the code versions:

for (( i=0 ; i < ${#VERSIONS[@]} ; ++i )); do

    version=${VERSIONS[$i]}
    options=${OPTIONS[$i]}

    echo
    echo
    echo "# Code Version: $version"
    echo "# ------------"
    echo

    # Create install directory
    install_dir="$install_prefix/$version"
    run "mkdir -p $install_dir" || die "Could not create $install_dir"

    # Loop over the packages.
    for package in ${PACKAGES[@]};  do
        echo
        echo "#    Package: $package"
        echo "#    -------"
        echo
        
        source_dir="$source_prefix/source/$package"
        build_dir="$build_prefix/$version/${package:0:1}"
        if test -d ${build_dir}; then
            run "rm -rf ${build_dir}"
        fi

        run "mkdir -p $build_dir" || die "Could not create directory $build_dir."
        run "cd $build_dir"
        echo "CMAKE_SYSTEM_NAME:STRING=Catamount" > $build_dir/CMakeCache.txt
        echo "DRACO_LIBRARY_TYPE:STRING=STATIC" >> $build_dir/CMakeCache.txt
        echo "CMAKE_C_COMPILER:FILEPATH=cc" >> $build_dir/CMakeCache.txt
        echo "CMAKE_CXX_COMPILER:FILEPATH=CC " >> $build_dir/CMakeCache.txt
        echo "CMAKE_Fortran_COMPILER:FILEPATH=ftn" >> $build_dir/CMakeCache.txt
        echo "MPIEXEC:FILEPATH=/usr/bin/aprun" >> $build_dir/CMakeCache.txt
        echo "MPIEXEC_NUMPROC_FLAG:STRING=-n" >> $build_dir/CMakeCache.txt
        echo "MPI_C_LIBRARIES:FILEPATH=" >> $build_dir/CMakeCache.txt
        echo "MPI_CXX_LIBRARIES:FILEPATH=" >> $build_dir/CMakeCache.txt
        echo "MPI_Fortran_LIBRARIES:FILEPATH=" >> $build_dir/CMakeCache.txt
        echo "MPI_C_INCLUDE_PATH:PATH=" >> $build_dir/CMakeCache.txt
        echo "MPI_CXX_INCLUDE_PATH:PATH=" >> $build_dir/CMakeCache.txt
        echo "MPI_Fortran_INCLUDE_PATH:PATH=" >> $build_dir/CMakeCache.txt
        echo "CMAKE_INSTALL_PREFIX:PATH=$install_dir" >> $build_dir/CMakeCache.txt
        sleep 2
        run "cat $build_dir/CMakeCache.txt"
        run "cmake \
            $options $CONFIG_BASE $source_dir" \
            || die "Could not configure in $build_dir from source at $source_dir"
        run "$MAKE all"  || die "Could not build code/tests in $build_dir"
        # only run tests for non-nr versions
#        case $version in
#        *_nr) ;;
#        *)    run "ctest" ;;
#        esac
        run "$MAKE install"  || die "Could not build code/tests in $build_dir"
        run "chmod -R $build_permissions $build_dir"

    done

    # Set access to install dir.
    run "chmod -R $install_permissions $install_dir"

done

# Set access to top level install dir.
run "chmod $install_permissions $install_prefix"
run "chgrp -R othello $install_prefix"

