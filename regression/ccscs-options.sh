#!/bin/bash
## -*- Mode: sh -*-
##---------------------------------------------------------------------------##
## File  : regression/sripts/ccscs-options.sh
## Date  : Tuesday, May 16, 2017, 21:42 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##
##
## Summary: Capture ccscs-specific features that plug into regression-master.sh.
##---------------------------------------------------------------------------##

# Main Options
export machine_name_long="Linux64 on CCS LAN"
platform_extra_params="clang coverage fulldiagnostics gcc630 nr perfbench valgrind"
pem_match=`echo $platform_extra_params | sed -e 's/[ ]/|/g'`

##---------------------------------------------------------------------------##
## End ccscs-options.sh
##---------------------------------------------------------------------------##
