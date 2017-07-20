#!/bin/bash

## -*- Mode: sh -*-
##---------------------------------------------------------------------------##
## File  : regression/metrics_report.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##
## Generate a code-coverage and LOC report and send it by email.
##---------------------------------------------------------------------------##

# Collect montly metrics for Jayenne codes and email a report.
# Use:
#    <path>/metrics_report.sh -e "email1 [email2]" -p "project1 [project2]"

# Typical use:
# ./metrics_report.sh -e kgt@lanl.gov -p "draco jayenne capsaicin"
# ./metrics_report.sh -e "jsbrock@lanl.gov jomc@lanl.gov gshipman@lanl.gov draco@lanl.gov" -p "draco"

##---------------------------------------------------------------------------##
## Support functions
##---------------------------------------------------------------------------##
print_use()
{
    echo " "
    echo "Combine code metric data and email report to recipients."
    echo " "
    echo "Usage: $0 -e \"email1 [email2]\" -p \"project1 [project2]\""
    echo " "
    echo "   [projects] = { draco, jayenne, capsaicin }."
    echo "   - Multiple projects or emails must be space delimeted and"
    echo "     in quotes."
}

##---------------------------------------------------------------------------##
## Sanity Checks and Setup
##---------------------------------------------------------------------------##

if [[ ! ${4} ]]; then
   echo "ERROR: You must provide at least 4 arguments."
   print_use
   exit 1
fi

mach=`uname -n`
if test "$mach" != "ccscs7.lanl.gov"; then
   echo "FATAL ERROR: This script must be run from ccscs7 (or the machine that"
   echo "             has the regression log files)."
   exit 1
fi

if [[ ! ${USER} ]]; then
   echo "FATAL ERROR: ENV{USER} not set.  Contact Kelly Thompson <kgt@lanl.gov>."
   exit 1
fi

# Store log and temporary files here.
logdir=/scratch/${USER}
work_dir=/scratch/regress
CLOC=/scratch/vendors/bin/cloc

# Ensure the work directory exists
if ! test -d $logdir; then
   mkdir -p $logdir || exit 1
fi

# Remove any exising log file.
logfile=$logdir/metrics.log
if test -f $logfile; then
   rm $logfile
fi
touch $logfile

if ! test -x $CLOC; then
   echo "FATAL ERROR: 'cloc' not found in PATH. Contact Kelly Thompson <kgt@lanl.gov>."
   exit 1
fi


if ! test -d $work_dir/cdash; then
   echo "FATAL ERROR: Regression work directory (work_dir=$work_dir/cdash) not found.  Contact Kelly Thompson <kgt@lanl.gov>."
   exit 1
fi


##---------------------------------------------------------------------------##
# Environment setup
##---------------------------------------------------------------------------##

# load some common bash functions
export rscriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [[ -f $rscriptdir/scripts/common.sh ]]; then
  source $rscriptdir/scripts/common.sh
else
  echo " "
  echo "FATAL ERROR: Unable to locate Draco's bash functions: "
  echo "   looking for .../regression/scripts/common.sh"
  echo "   searched rscriptdir = $rscriptdir"
  exit 1
fi



# Modules
# ----------------------------------------

if [[ `fn_exists module` == 1 ]]; then
  echo " "
  if [[ `declare -f module | grep -c LMOD` == 0 ]]; then
    # we have tcl modules
    echo "Found Tcl modules:"
    run "module list"
    echo "unloading"
    run "module purge"
    run "module list"
    unset dracomodules
    unset NoModules
    unset _LMFILES_
    unset MODULEPATH
    unset LOADEDMODULES
    unset MODULESHOME

    # If TCL modulefiles, switch to Lmod
    echo " "
    echo "Switching to Lmod modulefiles..."
    export MODULE_HOME=/scratch/vendors/spack.20170502/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/lmod-7.4.8-oytncsoih2sa4jdogz2ojvwly6mwle4n
    source $MODULE_HOME/lmod/lmod/init/bash || die "Can't find /mod/init/bash"
    run "module use /scratch/vendors/spack.20170502/share/spack/lmod/linux-rhel7-x86_64/Core" || die "Can't find Lmod Core modulefiles."

  else
    # we have Lmod modules
    echo "Found Lmod modules:"
    run "module avail"
    run "module list"
    run "module purge"
    run "module list"
  fi
else
  echo " "
  echo "No modules available"
  echo "Loading Lmod modulefiles..."
  export MODULE_HOME=/scratch/vendors/spack.20170502/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/lmod-7.4.8-oytncsoih2sa4jdogz2ojvwly6mwle4n
  source $MODULE_HOME/lmod/lmod/init/bash || die "Can't find /mod/init/bash"
  run "module use /scratch/vendors/spack.20170502/share/spack/lmod/linux-rhel7-x86_64/Core" || die "Can't find Lmod Core modulefiles."
fi

# Establish environment via Lmod...

# eospac, ndi, csk
run "module use --append /scratch/vendors/Modules.lmod"
run "module load bullseyecoverage"

COVDIR=`which covdir`
if ! test -x $COVDIR; then
   echo "FATAL ERROR: covdir not found in PATH  Contact Kelly Thompson <kgt@lanl.gov>."
   exit 1
fi

##---------------------------------------------------------------------------##
# Process arguments
##---------------------------------------------------------------------------##

while getopts ":e:p:" opt; do
case $opt in
e)
   recipients=$OPTARG
   #echo "This report will be emailed to $recipients"
   ;;
p)
   projects=$OPTARG
   #echo "Processing project(s): $projects"
   ;;
:)
   echo "Option -$OPTARG requires an argument."
   print_use
   exit 1
   ;;
\?)
   print_use
   exit 1
   ;;
esac
done

# redirect all script output to file
echo "The metric report will be saved to $logfile "
echo "and emailed to ${recipients}."
exec 1>${logfile}
exec 2>&1

##---------------------------------------------------------------------------##
# Banner
##---------------------------------------------------------------------------##

# defaults:
project_name="${projects}"
pp=""
ntrim_lines=7 # number of lines to trim from the end of cloc output.
# special cases:
if test "${projects}" = "draco"; then
   project_name="Draco"
   let ntrim_lines+=1
elif test "${projects}" = "draco jayenne"; then
   project_name="Jayenne"
   pp=" - Combined report for Draco and Jayenne"
   let ntrim_lines+=2
elif test "${projects}" = "draco capsaicin"; then
   project_name="Capsaicin"
   pp=" - Combined report for Draco and Capsaicin"
   let ntrim_lines+=2
elif test "${projects}" = "draco jayenne capsaicin"; then
   project_name="Jayenne and Capsaicin"
   pp=" - Combined report for Draco, Jayenne and Capsaicin"
   let ntrim_lines+=3
fi
# we need a negative number
ntrim_lines=`expr $ntrim_lines \* -1`

subj=`echo -n "Code Metrics for ${project_name}, "; date`

echo "======================================================================"
echo "${subj}"
echo "======================================================================"
echo " "
echo " "
echo "--------------------"
echo "${project_name} $pp"
echo "--------------------"
echo " "
# Lines of code report for $projects

echo "Lines of code"
echo "-------------"
cmd="${CLOC} --sum-reports --force-lang-def=/scratch/regress/draco/regression/cloc-lang.defs"
for proj in $projects; do
   if test -f ${work_dir}/cdash/${proj}/Nightly_gcc/Coverage/build/lines-of-code.log; then
     cmd="$cmd ${work_dir}/cdash/${proj}/Nightly_gcc/Coverage/build/lines-of-code.log"
   else
     cmd="$cmd ${work_dir}/cdash/${proj}/Nightly_gcc-develop/Coverage/build/lines-of-code.log"
   fi
done
# Use grep and head to clean up the output:
cmd="$cmd | grep -v sourceforge | grep -v AlDanial | head -n $ntrim_lines"
eval $cmd

echo " "
echo "Code coverage"
echo "-------------"
echo " "
export COVFILE=$logdir/metrics_report.cov
if test -f $COVFILE; then
   rm -f $COVFILE
fi
cmd="covmerge -q --mp --no-banner -c -f $COVFILE "
for proj in $projects; do
   if test -f ${work_dir}/cdash/${proj}/Nightly_gcc/Coverage/build/CMake.cov.bak; then
     cmd="$cmd ${work_dir}/cdash/${proj}/Nightly_gcc/Coverage/build/CMake.cov.bak "
   else
     cmd="$cmd ${work_dir}/cdash/${proj}/Nightly_gcc-develop/Coverage/build/CMake.cov.bak "
   fi
done
# create the new coverage file via covmerge
eval $cmd
# run covdir to generate a report (but omit entry for /source/src/)
olddir=`pwd`
cd $work_dir
export COVDIRCFG=$work_dir/draco/regression/mcovdir.cfg

# 1. Generate a merged report, then
# 2. Trim leading directory names.
# 3. Strip out the line '../source/src/'
# 4. Drop lines that are package subdirectories since the function coverate
#    numbers are already included in the parent.
covdir -w120 \
| sed -r 's%../source/src/([A-Za-z0-9+_/]+)/%\1/             %' \
| grep -v "../source/src/" \
| grep -v '[a-z]/[a-z]'

#| sed -e 's/Directory               /Directory/' \
#| sed -e 's/------------------------------------------------/---------------------------------/' \
#| sed -e 's/Total               /Total/' \
#| grep -v "source/.*src/ " | grep -v "source/ " | grep -v "^src"

echo " "
echo "* C/D Coverage is condition/decision coverage"
echo "  http://www.bullseye.com/coverage.html#basic_conditionDecision"

cd $olddir

# Send the email
/bin/mailx -r "${USER}@lanl.gov" -s "${subj}" ${recipients} < ${logfile}
