#!/bin/bash -l

if [[ `uname -n` =~ "darwin-fe" ]]; then
  echo "FATAL ERROR: This script must be run from a back-end node"
  exit 1
fi

# export dry_run=1
export steps="config build test"
echo -e "\nConfigure, Build and Test $buildflavor-$rttversion version of $package.\n"
$draco_script_dir/release.msub

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
