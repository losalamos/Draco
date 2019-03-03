#!/bin/bash
## -*- Mode: sh -*-
##---------------------------------------------------------------------------##
## File  : regression/sripts/cts1-options.sh
## Date  : Tuesday, May 16, 2017, 21:42 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2017, Triad National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##
##
## Summary: Capture cts1-specific features that plug into regression-master.sh.
##---------------------------------------------------------------------------##

# Main Options
machine=`uname -n`
export machine_class="cts1"
case $machine in
ba*)
    export machine_name_short="ba"
    export machine_name_long="Badger" ;;
gr*)
    export machine_name_short="gr"
    export machine_name_long="Grizzly" ;;
ic*)
    export machine_name_short="ic"
    export machine_name_long="Ice" ;;
fi*)
    export machine_name_short="fi"
    export machine_name_long="Fire" ;;
sn*)
    export machine_name_short="sn"
    export machine_name_long="Snow" ;;
esac
platform_extra_params="fulldiagnostics gcc610 gcc640 newtools nr perfbench valgrind vtest"
pem_match=`echo $platform_extra_params | sed -e 's/[ ]/|/g'`

##---------------------------------------------------------------------------##
## End cts1-options.sh
##---------------------------------------------------------------------------##
