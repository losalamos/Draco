#!/bin/bash

# Collect montly metrics for Jayenne codes and email a report.
# Use:
#    <path>/metrics_report.sh -e "email1 [email2]" -p "project1 [project2]"

# Typical use:
# ./metrics_report.sh -e kellyt@lanl.gov -p "draco clubimc wedgehog"
# ./metrics_report.sh -e "kgt@lanl.gov jhchang@lanl.gov" -p "draco capsaicin"
# ./metrics_report.sh -e "jsbrock@lanl.gov barcher@lanl.gov caw@lanl.gov jayenne@lanl.gov" -p "draco clubimc wedgehog"

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
    echo "   [projects] = { draco, clubimc, wedgehog, milagro, capsaicin }."
    echo "   - Multiple projects or emails must be space delimeted and"
    echo "     in quotes."
}

##---------------------------------------------------------------------------##
## Sanity Checks and Setup
##---------------------------------------------------------------------------##

if test "${4}x" = "x"; then
   echo "ERROR: You must provide at least 4 arguments."
   print_use
   exit 1
fi

mach=`uname -n`
if test "$mach" != "ccscs8.lanl.gov"; then
   echo "FATAL ERROR: This script must be run from ccscs8."
   exit 1
fi

if test "${USER}x" == "x"; then
   echo "FATAL ERROR: ENV{USER} not set.  Contact Kelly Thompson <kgt@lanl.gov>."
   exit 1
fi

# Ensure the work directory exists
if ! test -d /var/tmp/${USER}; then
   mkdir -p /var/tmp/${USER}
fi

# Remove any exising log file.
logfile=/var/tmp/${USER}/metrics.log
if test -f $logfile; then
   rm $logfile
fi
if ! test -d /var/tmp/${USER}; then
    mkdir /var/tmp/${USER}
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

CLOC=/home/regress/draco/regression/cloc
work_dir=/var/tmp/regress

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
elif test "${projects}" = "draco clubimc wedgehog"; then
   project_name="Jayenne"
   pp=" - Combined report for Draco, ClubIMC and Wedgehog"
elif test "${projects}" = "draco clubimc wedgehog milagro"; then
   project_name="Jayenne (with Milagro)"
   pp=" - Combined report for Draco, ClubIMC, Wedgehog and Milagro"
elif test "${projects}" = "draco capsaicin"; then
   project_name="Capsaicin"
   pp=" - Combined report for Draco and Capsaicin"
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
cmd="${CLOC} --sum-reports --read-lang-def=/home/regress/draco/regression/cloc-lang.defs"
for proj in $projects; do
   cmd="$cmd ${work_dir}/${proj}/Nightly_gcc/Coverage/build/lines-of-code.log "
done
# Use grep and head to clean up the output:
cmd="$cmd | grep -v sourceforge | head -n 24"
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
export COVFILE=`pwd`/metrics_report.cov
export COVDIRCFG=/home/regress/draco/regression/covdir.cfg
if test -f $COVFILE; then
   rm -f $COVFILE
fi
cmd="covmerge -q --mp --no-banner -c -f $COVFILE "
for proj in $projects; do
   cmd="$cmd ${work_dir}/${proj}/Nightly_gcc/Coverage/build/CMake.cov "
done
# create the new coverage file via covmerge
eval $cmd
# run covdir to generate a report
covdir

echo " "
echo "* C/D Coverage is condition/decision coverage"
echo "  http://www.bullseye.com/coverage.html#basic_conditionDecision"

# Send the email

#echo " "
#echo /bin/mailx -s \"${subj}\" ${recipients} < ${logfile}
/bin/mailx -s "${subj}" ${recipients} < ${logfile}
#/bin/mailx -s "${subj}" kellyt@lanl.gov < ${logfile}
#/bin/mailx -s "${subj}" jsbrock@lanl.gov barcher@lanl.gov jayenne@lanl.gov < ${logfile}
