#! /usr/bin/sh
#  ========================================================================  #
# 
#       Author: Randy Roberts 
#               Los Alamos National Laboratory
#       Date:   Fri Oct  8 17:20:38 MDT 1999
# 
#       Copyright (c) 1999 U. S. Department of Energy. All rights reserved.
# 
#       $Id$
# 
#  NAME
#       get_include_tree.sh - Finds all the include files that are
#                             used inside a package
#  SYNOPSIS
#       get_include_tree.sh  dir
#  DESCRIPTION
#       get_include_tree.sh looks through the .d dependency files generated
#       as a result of a configure and make, finding the package and
#       files names that are used.
#
#  ========================================================================  #

# Look through the depend files and extract the .hh files
dependFiles=`find $1 -name '*.d' -print`
hhFiles=`cat $dependFiles | grep '.hh' | awk '{print $4}' | sort | uniq`

# Loop though the .hh files to get the package and file name
for file in $hhFiles
do
    if echo $file| grep include >/dev/null
    then
       package=`echo $file | sed -n '/include/s#.*/include/\([^/]*\)/.*#\1#p'`
       name=`basename $file '.hh'`
       echo "$package/$name"
       # echo $file
    fi
done |
sort | uniq

