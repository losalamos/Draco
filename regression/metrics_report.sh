#!/bin/bash

# Collect montly metrics for Jayenne codes and email a report.

# Remove any exising log file.
logfile=/var/tmp/kellyt/metrics.log
if test -f $logfile; then
   rm $logfile
fi

# Environment setup

if test -z "$MODULESHOME"; then
  # This is a new login
  if test -f /home/regress/environment/Modules/init/bash; then
    source /home/regress/environment/Modules/init/bash
  fi
fi
module load bullseyecoverage/8.4.12

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
cd /home/regress/cmake_draco/Nightly_gcc/Coverage/build
export COVFILE=`pwd`/CMake.cov
export COVDIRCFG=`pwd`/covclass_cmake.cfg
cd ../source/src
covdir

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
export COVFILE=`pwd`/CMake.cov
export COVDIRCFG=`pwd`/covclass_cmake.cfg
cd ../source/src
covdir

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
export COVFILE=`pwd`/CMake.cov
export COVDIRCFG=`pwd`/covclass_cmake.cfg
cd ../source/src
covdir

# Jayenne - Milagro

# echo " "
# echo " "
# echo "--------------------"
# echo "Jayenne - Milagro"
# echo "--------------------"
# echo " "
# echo "Lines of code"
# echo "-------------"
# cat /home/regress/cmake_jayenne/milagro/Nightly_gcc/Coverage/build/lines-of-code.log
# echo " "
# echo "Code coverage"
# echo "-------------"
# cd /home/regress/cmake_jayenne/milagro/Nightly_gcc/Coverage/build
# export COVFILE=`pwd`/CMake.cov
# export COVDIRCFG=`pwd`/covclass_cmake.cfg
# cd ../source/src
# covdir

 # svn diff -r"{2012-01-09}":"{2012-07-09}" --summarize &> ../m_mod_files.log
 #  num_files_changed=`cat ../m_mod_files.log | awk '{print $2}' | sort -u | wc -l`
 #   echo $num_files_changed 
 #  svn log -r"{2012-01-09}":"{2012-07-09}" --verbose &> ../m_msg.log
 #  number_of_commits=`cat ../m_msg.log | grep "Changed paths:" | wc -l`
 #  echo $number_of_commits
 #  files=`find . -name '*.hh' -o -name '*.cc' -o -name '*.txt' -o -name '*.cmake' -o -name '*.in' -o -name '*.h'`
 #  svn annotate $files > ../m_file_list
 #  user_list=`cat ../m_file_list | awk '{print $2}' | sort -u`
 #   echo $user_list 
 #  for name in $user_list; do numlines=`grep $name ../m_file_list | wc -l`; echo "$numlines: $name"; done > ../author_loc
 #   cat ../author_loc | sort -rn

# Send the email

/bin/mailx -s "${subj}" kellyt@lanl.gov < ${logfile}
#/bin/mailx -s "${subj}" jsbrock@lanl.gov barcher@lanl.gov jayenne@lanl.gov < ${logfile}

# jayenne@lanl.gov
# jsbrock@lanl.gov
# barcher@lanl.gov
