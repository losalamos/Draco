#!/bin/bash
## -*- Mode: sh -*-
##---------------------------------------------------------------------------##
## File  : regression/sripts/ml-options.sh
## Date  : Tuesday, May 16, 2017, 21:42 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##
##
## Summary: Capture ml-specific features that plug into regression-master.sh.
##---------------------------------------------------------------------------##

# Main Options
machine_name_long="Moonlight"
platform_extra_params="fulldiagnostics nr perfbench pgi valgrind"
pem_match=`echo $platform_extra_params | sed -e 's/[ ]/|/g'`

##---------------------------------------------------------------------------##
## End ml-options.sh
##---------------------------------------------------------------------------##
