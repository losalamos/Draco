#!/bin/bash

# Collect montly metrics for Jayenne codes and email a report.
# Use:
#    <path>/metrics_report.sh -e "email1 [email2]" -p "project1 [project2]"

# Typical use:
# ./metrics_report.sh -e kgt@lanl.gov -p "draco jayenne capsaicin"
# ./metrics_report.sh -e "jsbrock@lanl.gov jomc@lanl.gov sriram@lanl.gov draco@lanl.gov" -p "draco jayenne capsaicin"

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
   echo "FATAL ERROR: This script must be run from ccscs7."
   exit 1
fi

if [[ ! ${USER} ]]; then
   echo "FATAL ERROR: ENV{USER} not set.  Contact Kelly Thompson <kgt@lanl.gov>."
   exit 1
fi

# Ensure the work directory exists
if ! test -d /scratch/${USER}; then
   mkdir -p /scratch/${USER}
fi

# Remove any exising log file.
logfile=/scratch/${USER}/metrics.log
if test -f $logfile; then
   rm $logfile
fi
touch $logfile

##---------------------------------------------------------------------------##
# Environment setup
##---------------------------------------------------------------------------##

if test -z "$MODULESHOME"; then
  # This is a new login
  if test -f /ccs/codes/radtran/vendors/modules-3.2.9/init/bash; then
    source /ccs/codes/radtran/vendors/modules-3.2.9/init/bash
  fi
fi
module load bullseyecoverage

CLOC=/scratch/vendors/bin/cloc
work_dir=/scratch/regress/cdash

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
# special cases:
if test "${projects}" = "draco"; then
   project_name="Draco"
elif test "${projects}" = "draco jayenne"; then
   project_name="Jayenne"
   pp=" - Combined report for Draco and Jayenne"
elif test "${projects}" = "draco capsaicin"; then
   project_name="Capsaicin"
   pp=" - Combined report for Draco and Capsaicin"
elif test "${projects}" = "draco jayenne capsaicin"; then
   project_name="Jayenne and Capsaicin"
   pp=" - Combined report for Draco, Jayenne and Capsaicin"
fi

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
   if test -f ${work_dir}/${proj}/Nightly_gcc/Coverage/build/lines-of-code.log; then
     cmd="$cmd ${work_dir}/${proj}/Nightly_gcc/Coverage/build/lines-of-code.log"
   else
     cmd="$cmd ${work_dir}/${proj}/Nightly_gcc-develop/Coverage/build/lines-of-code.log"
   fi
done
# Use grep and head to clean up the output:
cmd="$cmd | grep -v sourceforge | head -n 28 | grep -v AlDanial"
eval $cmd

# ${CLOC} --sum-reports \
# ${work_dir}/draco/Nightly_gcc/Coverage/build/lines-of-code.log \
# ${work_dir}/clubimc/Nightly_gcc/Coverage/build/lines-of-code.log \
# ${work_dir}/wedgehog/Nightly_gcc/Coverage/build/lines-of-code.log \
# | grep -v sourceforge | head -n 24

echo " "
echo "Code coverage"
echo "-------------"
echo " "
export COVFILE=/scratch/$USER/metrics_report.cov
if test -f $COVFILE; then
   rm -f $COVFILE
fi
cmd="covmerge -q --mp --no-banner -c -f $COVFILE "
for proj in $projects; do
   if test -f ${work_dir}/${proj}/Nightly_gcc/Coverage/build/CMake.cov.bak; then
     cmd="$cmd ${work_dir}/${proj}/Nightly_gcc/Coverage/build/CMake.cov.bak "
   else
     cmd="$cmd ${work_dir}/${proj}/Nightly_gcc-develop/Coverage/build/CMake.cov.bak "
   fi
done
# create the new coverage file via covmerge
eval $cmd
# run covdir to generate a report (but omit entry for /source/src/)
olddir=`pwd`
cd $work_dir/..
export COVDIRCFG=/scratch/regress/draco/regression/mcovdir.cfg
covdir -w120 | sed -r \
's%../../../scratch/source/src/([A-Za-z0-9+_]+)/%\1/             %' \
| sed -e 's/Directory               /Directory/' \
| sed -e 's/------------------------------------------------/---------------------------------/' \
| sed -e 's/Total               /Total/' \
| grep -v "../source/src/"
#| grep -v "source/.*src/ " | grep -v "source/ " | grep -v "^src"

echo " "
echo "* C/D Coverage is condition/decision coverage"
echo "  http://www.bullseye.com/coverage.html#basic_conditionDecision"

cd $olddir

# Send the email
/bin/mailx -r "${USER}@lanl.gov" -s "${subj}" ${recipients} < ${logfile}
