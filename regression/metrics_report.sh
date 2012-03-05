#!/bin/bash

# Collect montly metrics for Jayenne codes and email a report.

# Remove any exising log file.
logfile=/var/tmp/kellyt/metrics.log
if test -f $logfile; then
   rm $logfile
fi

# redirect all script output to file
exec 1>${logfile}
exec 2>&1

# Banner

subj=`echo -n "Code Metrics for Draco and Jayenne, "; date`

echo "======================================================================"
echo ${subj}
echo "======================================================================"

# Draco

echo " "
echo "--------------------"
echo "Draco"
echo "--------------------"
echo " "
echo "Lines of code"
echo "-------------"
cat /home/regress/cmake_draco/Nightly_gcc/Coverage/build/lines-of-code.log
echo " "
echo "Code coverage"
echo "-------------"
module load bullseyecoverage/8.4.12
cd /home/regress/cmake_draco/Nightly_gcc/Coverage/build
COVFILE=`pwd`/CMake.cov
COVDIRCFG=`pwd`/covclass_cmake.cfg
cd ../source/src; covdir

# Jayenne - ClubIMC

echo " "
echo " "
echo "--------------------"
echo "Jayenne - ClubIMC"
echo "--------------------"
echo " "
echo "Lines of code"
echo "-------------"
cat /home/regress/cmake_jayenne/clubimc/Nightly_gcc/Coverage/build/lines-of-code.log
echo " "
echo "Code coverage"
echo "-------------"
cd /home/regress/cmake_jayenne/clubimc/Nightly_gcc/Coverage/build
COVFILE=`pwd`/CMake.cov
COVDIRCFG=`pwd`/covclass_cmake.cfg
cd ../source/src; covdir

# Jayenne - Wedgehog

echo " "
echo " "
echo "--------------------"
echo "Jayenne - Wedgehog"
echo "--------------------"
echo " "
echo "Lines of code"
echo "-------------"
cat /home/regress/cmake_jayenne/wedgehog/Nightly_gcc/Coverage/build/lines-of-code.log
echo " "
echo "Code coverage"
echo "-------------"
cd /home/regress/cmake_jayenne/wedgehog/Nightly_gcc/Coverage/build
COVFILE=`pwd`/CMake.cov
COVDIRCFG=`pwd`/covclass_cmake.cfg
cd ../source/src; covdir

# Send the email

/bin/mailx -s "${subj}" jayenne@lanl.gov < ${logfile}

