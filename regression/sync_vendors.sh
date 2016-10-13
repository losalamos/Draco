#!/bin/bash

# Run from ccscs7:
#   - Mirror ccscs7:/scratch/vendors -> /ccs/codes/radtran/vendors/rhel72
#     Do not mirror the ndi directory (300+ GB)
#   - Mirror ccscs7:/scratch/vendors  -> ccscs[234568]:/scratch/vendors

function run ()
{
    echo $1;
    if ! [ $dry_run ]; then
        eval $1;
    fi
}

machine=`uname -n`
vdir=/scratch/vendors
r72v=/ccs/codes/radtran/vendors/rhel72vendors
ccs_servers="2 3 4 5 6 7 8"

echo "Rsync vendor directory to /ccs/codes/radtran/vendors and to"
echo "/scratch/vendors on:"
for m in $ccs_servers; do echo " - ccscs${m}"; done

# make a backup copy of vendors to $r72v

if ! test -d $vdir; then
  echo "Source directory $vdir is missing."
  exit 1
fi

#cd $vdir
#find . -name '*~' -exec echo rm -f {} \;
echo " "
echo "Clean up permissions on source files..."
run "chgrp -R draco $vdir"
run "chmod -R g+rwX,o=g-w $vdir"

echo " "
echo "Save a copy on /ccs/codes/radtran..."
run "rsync -av --exclude 'ndi' --delete $vdir/ $r72v"
# -vaum

# rsync vendors ccscs7 -> other machines.
# limit network to 50 MB/sec (400 mbps)
echo " "
echo "Rsync $vdir to other ccs-net servers... "
for m in $ccs_servers; do
  if test `uname -n | grep $m | wc -l` = 0; then
    case $m in
      5)
        # do not copy ndi, not enough space
        run "rsync -av --exclude ndi --delete --bwlimit=50000 $vdir/ ccscs${m}:$vdir"
        ;;
      *)
        run "rsync -av --delete --bwlimit=50000 $vdir/ ccscs${m}:$vdir"
        # -vaum
        ;;
    esac
  fi
done

echo " "
echo "done"

# Rsync directories
# cd $mastervdir
# dirs=`\ls -1 $mastervdir`
# for dir in $dirs; do
#     shortdir=`echo $dir | sed -e 's/-.*//'`

#     case $shortdir in
#     modules* | Modules* | deprecated | environment | sources | win32 )
#        # do not process
#        ;;
#     emacs* )
#        # do not process
#        ;;
#     *)
#        echo " "
#        echo "cd $masterdir"
#        echo "rsync -vaum $dir ${VENDOR_DIR}"
#        rsync -vaum $dir ${VENDOR_DIR}
#        ;;
#     esac
# done

# chgrp -R draco ${VENDOR_DIR}
# chmod -R g+rwX,o=g-w ${VENDOR_DIR}
