#! /bin/sh

# This script descends subdirectories looking for Release.cc.
# It then compares the pkg_release string to existing CVS tags.

# Usage: releaseCheck.sh

doit()
{
    file=$1
    release=`sed -n -e "/@(#)/ s/.*@(#)\(.*\)\";/\1/p" < $file`
    if [ "$release" = "" ]
    then
       echo $file pkg_release is ill-formed or non-existant
       return 1
    fi
    workingRevision=`cvs status -v $file 2>/dev/null |
                     egrep "Working revision:" | awk '{print $3}' | sed 's/)//'`

    if cvs status -v $file 2>/dev/null | egrep last_stable >/dev/null
    then
       last_stable=`cvs status -v $file 2>/dev/null |
                     egrep "last_stable" | awk '{print $3}' | sed 's/)//'`
    fi

    if cvs status -v $file 2>/dev/null | egrep $release >/dev/null
    then
       releaseRev=`cvs status -v $file 2>/dev/null |
                     egrep "$release" | awk '{print $3}' | sed 's/)//'`
       echo $file: release=$release tag found
       echo "\t" release: $releaseRev working: $workingRevision last_stable: $last_stable
    else
       echo $file: release=$release tag not found
    fi
}

files=`find . -name Release.cc -print`

for f in $files
do
  doit $f
done

