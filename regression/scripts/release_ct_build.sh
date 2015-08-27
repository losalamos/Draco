#!/bin/bash -l

#----------------------------------------------------------------------#
# The script starts here
#----------------------------------------------------------------------#

# Helpful functions:
die () { echo "ERROR: $1"; exit 1;}

run () {
    echo $1
    if ! test $dry_run; then
    eval $1
    fi
}

# Permissions - new files should be marked u+rwx,g+rwx,o+rx
umask 0002
if test `groups | grep othello | wc -l` = 1; then
    install_group="othello"
    install_permissions="g+rwX,o-rwX"
else
    install_group="draco"
    install_permissions="g+rwX,o=g-w"
fi
build_permissions="g+rwX"

# environment (use draco modules)
run "module use /usr/projects/draco/vendors/Modules/hpc"
run "module load friendly-testing user_contrib"
run "module unload ndi ParMetis SuperLU_DIST trilinos"
run "module unload lapack gsl intel"
run "module unload cmake numdiff svn"
run "module unload PrgEnv-intel PrgEnv-pgi"
run "module unload papi perftools"
run "module load PrgEnv-intel"
run "module unload xt-libsci xt-totalview"
# run "module swap intel intel/14.0.4.211"
run "module swap intel intel/15.0.3"
run "module load gsl/1.15"
run "module load cmake/3.3.1 numdiff svn"
run "module load trilinos SuperLU_DIST"
run "module load ParMetis ndi random123 eospac/v6.2.4"
run "module list"

export OMP_NUM_THREADS=8
export SCRATCH="/lscratch1"
export MAKEOPTS="-j 48 -l 48"
export CTEST="-j8"

# Define your source and build information here.

ddir="draco-6_17_0"
platform="ct"
dmpi=craympich2 # mpt string?
df90=intel1503
dcpp=intel1503

source_prefix="/usr/projects/draco/$ddir"
build_prefix="$SCRATCH/$USER/$ddir/$platform-${dmpi}-${df90}-${dcpp}"
install_prefix="$source_prefix/$platform-${dmpi}-${df90}-${dcpp}"

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

        # echo "CMAKE_SYSTEM_NAME:STRING=Catamount" > $build_dir/CMakeCache.txt
        # echo "DRACO_LIBRARY_TYPE:STRING=STATIC" >> $build_dir/CMakeCache.txt
        # echo "CMAKE_C_COMPILER:FILEPATH=cc" >> $build_dir/CMakeCache.txt
        # echo "CMAKE_CXX_COMPILER:FILEPATH=CC " >> $build_dir/CMakeCache.txt
        # echo "CMAKE_Fortran_COMPILER:FILEPATH=ftn" >> $build_dir/CMakeCache.txt
        # echo "MPIEXEC:FILEPATH=/usr/bin/aprun" >> $build_dir/CMakeCache.txt
        # echo "MPIEXEC_NUMPROC_FLAG:STRING=-n" >> $build_dir/CMakeCache.txt
        # echo "MPI_C_LIBRARIES:FILEPATH=" >> $build_dir/CMakeCache.txt
        # echo "MPI_CXX_LIBRARIES:FILEPATH=" >> $build_dir/CMakeCache.txt
        # echo "MPI_Fortran_LIBRARIES:FILEPATH=" >> $build_dir/CMakeCache.txt
        # echo "MPI_C_INCLUDE_PATH:PATH=" >> $build_dir/CMakeCache.txt
        # echo "MPI_CXX_INCLUDE_PATH:PATH=" >> $build_dir/CMakeCache.txt
        # echo "MPI_Fortran_INCLUDE_PATH:PATH=" >> $build_dir/CMakeCache.txt
        # echo "CMAKE_INSTALL_PREFIX:PATH=$install_dir" >> $build_dir/CMakeCache.txt
        # sleep 2
        # run "cat $build_dir/CMakeCache.txt"
        run "cmake -C $source_dir/config/CrayConfig.cmake \
            $options $CONFIG_BASE -DCMAKE_INSTALL_PREFIX:PATH=$install_dir $source_dir" \
            || die "Could not configure in $build_dir from source at $source_dir"
        run "make $MAKEOPTS all"  || die "Could not build code/tests in $build_dir"
        # only run tests for non-nr versions
#        case $version in
#        *_nr) ;;
#        *)    run "$CTEST" ;;
#        esac
        run "make $MAKEOPTS install"  || die "Could not build code/tests in $build_dir"
        run "chmod -R $build_permissions $build_dir"

    done

    # Set access to install dir.
    run "chmod -R $install_permissions $install_dir"

done

# Set access to top level install dir.
run "chmod $install_permissions $install_prefix"
run "chgrp -R ${install_group} $install_prefix"
