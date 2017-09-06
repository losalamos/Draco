#!/bin/bash

# format_commit.sh <file>
file=$1
tooldir=`dirname $0`
tooldir=`(cd $tooldir; pwd)`

EMACS=/ccs/codes/radtran/vendors/emacs-24.4/bin/emacs
if test -x $EMACS; then
    cppfile=`echo $file | egrep -c -e "[.]hh|[.]cc"`
    if test "${cppfile}" = "1"; then
        # echo "$EMACS -batch $file -l $tooldir/format_commit.el -f save-buffer"
        $EMACS -batch $file -l $tooldir/format_commit.el -f save-buffer &> /dev/null
    fi
fi
