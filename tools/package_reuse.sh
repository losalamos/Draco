#! /usr/bin/sh
#  ========================================================================  #
# 
#       Author: Mark G. Gray, John M. McGhee
#               Los Alamos National Laboratory
#       Date:   Thu Oct 14 17:02:00 MDT 1999
# 
#       Copyright (c) 1999 U. S. Department of Energy. All rights reserved.
# 
#       $Id$
#
#  NAME
#        package_reuse.sh - Estimates the amount of code reuse that 
#        a group of packages is achieving
#                      
#  SYNOPSIS
#       package_reuse.sh [-s file] filelist....
#
#  DESCRIPTION
#       package_reuse.sh estimates the amount of code reuse that a
#       group of packages is achieving. "file" is a file containing
#       basenames and lines of code for the reuse library. If -s is not
#       present "file" defaults to "basename_sizes.txt". This file can
#       be produced with the utility "draco/tools/basename_loc".
#       "filelist" is a list of files containing include information
#       from systems using the reuse library. These files can be
#       produced with the utility "draco/tools/get_include_tree.sh".
#       Output consists of package name, total number of lines of
#       code used from the package by all the client systems, and 
#       the number of lines in the package. Reuse factor is estimated by
#       the ratio (total use)/(package size).
#       
#  ========================================================================  #

# Get the name of the file that contains lines of code per basename
if [ $1 = -s ]
then
size_file=$2 
shift 2
else
size_file="basename_sizes.txt"
fi

#Check to see if the size file exits
if [ ! -f $size_file ];
then
 echo $size_file " does not exist."
 exit 1
fi

#Produce the reuse estimate
echo "             Package" "   Total"   "  Actual"  " Reuse "
echo "                Name" "   Use  "   "  Size  "  " Factor"
echo "             -------" "   -----"   "  ------"  " ------"
cat $@ | awk '{tmp[$1]++ } END {for(i in tmp) print i, " ", tmp[i] }' | sort |
join $size_file - | awk '{gsub(/\/.*/,"",$1); tmp[$1]+=$2*$3; \
sum+=$2*$3; tmp2[$1]+=$2; sum2+=$2} END \
{for(i in tmp) { j=tmp[i]/tmp2[i]; k=10*(tmp[i]%tmp2[i]);  k=k/tmp2[i]; \
printf("%20s\t%d\t%d\t%d%1s%d\n", i, tmp[i], tmp2[i], j, ".", k); } \
printf("%20s\t%d\t%d\t%d%1s%d\n","Total", sum, sum2, sum/sum2, ".", \
(10*(sum%sum2))/sum2 )}' | sort -n +1

#  ==================== end of package_reuse.sh ===============================  #
