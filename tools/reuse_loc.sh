#! /usr/bin/sh
#  ========================================================================  #
# 
#       Author: Mark G. Gray
#               Los Alamos National Laboratory
#       Date:   Thu Oct 14 08:45:00 MDT 1999
# 
#       Copyright (c) 1999 U. S. Department of Energy. All rights reserved.
# 
#       $Id$
# 
#  NAME
#       reuse_loc - Given two input files produced by draco/tools/basename_loc
#                   and draco/tools/get_include_tree.sh, estimates the 
#                   lines of C++ code that a client system is using from a 
#                   library. 
#
#  SYNOPSIS
#       reuse_loc.sh  file1 file2
#
#  DESCRIPTION
#       draco/tools/basename_loc produces a listing of the lines of C++
#       code associated with each file basename in a reuse library
#       directory tree. 
#       draco/tools/get_include_tree.sh produces a listing of all the 
#       include files in a client system directory tree. 
#       reuse_loc.sh finds the intersection of these two data files, and
#       totals the lines of C++ code associated with the intersection.
#       This result is an estimate of how many lines of code a client
#       system is using from the reuse library.
#
#  ========================================================================  #

#Get command line arguments, if any
if [ $# != 2 ]
then
   echo "Error - Two filename arguments required"
   echo
   exit 1
fi

#Check for existence of the files to be scanned
if [ ! -f $1 ];
then
 echo "Error - file " $1 " does not exist."
 exit 1
fi

if [ ! -f $2 ];
then
 echo "Error - file " $2 " does not exist."
 exit 1
fi

#Scan the files to determine the intersection.
join  $1 $2 | awk '{sum+=$2} END {print sum}'


#  ====================== end of reuse_loc.sh =============================  #
