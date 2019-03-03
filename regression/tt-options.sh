#!/bin/bash
## -*- Mode: sh -*-
##---------------------------------------------------------------------------##
## File  : regression/sripts/tt-options.sh
## Date  : Tuesday, May 16, 2017, 21:42 pm
## Author: Kelly Thompson
## Note  : Copyright (C) 2017, Triad National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##
##
## Summary: Capture tt-specific features that plug into regression-master.sh.
##---------------------------------------------------------------------------##

# Main Options
export machine_class="tt"
export machine_name_short="tt"
export machine_name_long="Trinitite"
platform_extra_params="fulldiagnostics knl nr perfbench vtest"
pem_match=`echo $platform_extra_params | sed -e 's/[ ]/|/g'`

##---------------------------------------------------------------------------##
## End tt-options.sh
##---------------------------------------------------------------------------##
