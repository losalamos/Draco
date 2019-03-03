#!/bin/bash
## -*- Mode: sh -*-
##---------------------------------------------------------------------------##
## File  : regression/alist.sh
## Date  : Tuesday, May 31, 2016, 14:48 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2016-2019, Triad National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##
## Generate a list of authors by using git annotate
##---------------------------------------------------------------------------##

# Should be run from the draco top-level source directory.
if ! test -f CMakeLists.txt; then
  echo "ERROR: Must be run from the top level draco source directory."
  exit 1
fi
if ! test -f README.md; then
  echo "ERROR: Must be run from the top level draco source directory."
  exit 1
fi

run () {
   echo $1
   if ! test $dry_run; then
      eval $1
   fi
}

function math ()
{
    calc="${@//p/+}";
    calc="${calc//x/*}";
    bc -l <<< "scale=10;$calc"
}

#--------------------------------------------------------------------------------
# Generate basic statistics
# Exclude some files (especially binary and gold files).

# Generate a temporary file
author_loc=`mktemp`

# all files
tmp=`git ls-files`

# exclude data files
unset projectfiles
for file in $tmp; do
  case $file in
    *bench | *css | *dia | *eps | *fig | *ipcress | *jpg | *logfile )
      # drop these files from the count.
      ;;
    *mesh | *output.in | *png | *ps | *xs4 )
      # drop these files from the count.
      ;;
    *)
      if [[ -f $file ]]; then projectfiles="$projectfiles $file"; fi ;;
  esac
done

echo " "
echo "Collecting data.  This may take a few minutes..."
echo " "
echo $projectfiles | xargs -n1 git blame -w | perl -n -e '/^.*?\((.*?)\s+[\d]{4}/; print $1,"\n"' | sort -f | uniq -c | sort -rn > $author_loc

#--------------------------------------------------------------------------------
# Use pretty names and then de-duplicate.

# Put the user list into an array (some entries may have spaces).
IFS=$'\r\n' GLOBIGNORE='*' command eval 'entries=($(cat $author_loc))'

if [[ -f $author_loc ]]; then rm $author_loc; fi

for entry in "${entries[@]}"; do

  current_loc=`echo $entry | sed -e 's/[ ].*//'`
  current_author=`echo $entry | sed -e 's/[^ ]* //'`

  if [[ $current_loc == "" ]]; then current_loc=0; fi

  case $current_author in
    along | Alex*  ) current_author="Alex R. Long" ;;
    clevelam | cleveland | *Cleveland) current_author="Matt A. Cleveland" ;;
    gaber) current_author="Gabe M. Rockefeller" ;;
    hkpark) current_author="HyeongKae Park" ;;
    jdd) current_author="Jeff D. Densmore" ;;
    jhchang | Jae* ) current_author="Jae H. Chang" ;;
    keadyk | Kendra* ) current_author="Kendra P. Keady" ;;
    kellyt | kgt | 107638 | Kelly*) current_author="Kelly G. Thompson" ;;
    kgbudge) current_author="Kent G. Budge" ;;
    kwang) current_author="Katherine J. Wang" ;;
    lowrie | *Lowrie) current_author="Rob B. Lowrie" ;;
    lpritch) current_author="Lori A. Pritchett-Sheats" ;;
    maxrosa) current_author="Massimiliano Rosa" ;;
    ntmyers) current_author="Nick Meyers" ;;
    pahrens) current_author="Peter Ahrens" ;;
    talbotp) current_author="Paul Talbot" ;;
    tkelley | Tim*) current_author="Tim Kelley" ;;
    tmonster) current_author="Todd J. Urbatsch" ;;
    warsa | *Warsa) current_author="James S. Warsa" ;;
    wollaber) current_author="Allan B. Wollaber" ;;
    wollaege | *Wollaeger) current_author="Ryan T. Wollaeger" ;;
  esac

  echo "$current_author:$current_loc"

done | sort -rn > $author_loc

# cat $author_loc

#--------------------------------------------------------------------------------
# Merge data

prev_author="none"
prev_loc=0

# Put the user list into an array (some entries may have spaces).
IFS=$'\r\n' GLOBIGNORE='*' command eval 'entries=($(cat $author_loc))'
rm $author_loc
for entry in "${entries[@]}"; do

  current_author=`echo $entry | sed -e 's/:.*//'`
  current_loc=`echo $entry | sed -e 's/.*://'`

  if [[ $current_author == $prev_author ]]; then
    prev_loc=`math "$prev_loc + $current_loc"`
  else
    # previous author information before starting a new author.
    if ! [[ $prev_author == "none" ]]; then
      echo "${prev_loc}:${prev_author}" >> $author_loc
    fi

    prev_author=$current_author
    prev_loc=$current_loc
  fi

done
# the last entry
echo "${prev_loc}:${prev_author}" >> $author_loc

cat $author_loc | sort -rn

# cleanup
rm $author_loc $author_list_file

#------------------------------------------------------------------------------#
# end alist.sh
#------------------------------------------------------------------------------#
