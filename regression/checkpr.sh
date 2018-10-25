#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/checkpr.sh
## Date  : Wednesday, Mar 22, 2017, 16:01 pm
## Author: Kelly Thompson <kgt@lanl.gov>
## Note  : Copyright (C) 2017, Los Alamos National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

# Use:
#   <path>/checkpr.sh [options]
#
# Options:
#   -p <project>  - project name, {draco, jayenne, capasaicin}
#   -f <number>   - pr number
#   -h            - help message
#   -r            - special run mode that uses the regress account's
#                   credentials.
#   -t            - remove the last-draco tagfile (when was draco last built?).
#
# <number> must be an integer value that represents the pull request.

##---------------------------------------------------------------------------##
## Environment
##---------------------------------------------------------------------------##

# Enable job control
set -m

# load some common bash functions
export rscriptdir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" )
if ! [[ -d $rscriptdir ]]; then
  export rscriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
fi
if [[ -f $rscriptdir/scripts/common.sh ]]; then
  echo "source $rscriptdir/scripts/common.sh"
  source $rscriptdir/scripts/common.sh
else
  echo " "
  echo "FATAL ERROR: Unable to locate Draco's bash functions: "
  echo "   looking for .../regression/scripts/common.sh"
  echo "   searched rscriptdir = $rscriptdir"
  exit 1
fi

# defaults
export target=`uname -n | sed -e 's/[.].*//g'`
scratchdir=`selectscratchdir`
export logdir="$scratchdir/$USER/logs"
pr=develop
project=draco
regress_mode="off"

##---------------------------------------------------------------------------##
## Support functions
##---------------------------------------------------------------------------##

print_use()
{
  echo " "
  echo "Usage: ${0##*/} -h -p [draco|jayenne|capsaicin]"
  echo "       -f <git branch name> -r"
  echo " "
  echo "All arguments are optional,  The first value listed is the default value."
  echo "   -h    help           prints this message and exits."
  echo "   -f    git feature branch, default=develop"
  echo "         common: 'develop', '42'"
  echo "   -p    project name = { draco, jayenne, capsaicin }"
  echo "   -r    special run mode that uses the regress account's credentials."
  echo "   -t    remove the last-draco tagfile."
  echo " "
  echo "Examples:"
  echo "  ${0##*/} -p jayenne -f 42"
  echo "  ${0##*/} -p draco -t"
  echo "  ${0##*/} -p capsaicin -f develop"
  echo " "
}

##---------------------------------------------------------------------------##
## Command options
##---------------------------------------------------------------------------##

while getopts ":f:hp:rt" opt; do
case $opt in
f)  pr=$OPTARG ;;
h)  print_use; exit 0 ;;
p)  project=$OPTARG ;;
r)  regress_mode="on" ;;
t)  rmlastdracotag="on" ;;
\?) echo "" ;echo "invalid option: -$OPTARG"; print_use; exit 1 ;;
:)  echo "" ;echo "option -$OPTARG requires an argument."; print_use; exit 1 ;;
esac
done

##---------------------------------------------------------------------------##
## Sanity Checks for input
##---------------------------------------------------------------------------##

case $project in
  draco | jayenne | capsaicin ) # known projects, continue
    ;;
  *)  echo "" ;echo "FATAL ERROR: unknown project name (-p) = ${proj}"
    print_use; exit 1 ;;
esac

if [[ $pr =~ "pr" ]]; then
  echo -e "\nFATAL ERROR, the '-f' option expects a number (i.e.: no 'pr' prefix) or the string 'develop'."
  print_use;
  exit 1;
fi

if [[ $regress_mode == "on" ]]; then
  if ! [[ $LOGNAME == "kellyt" ]]; then
    echo ""; echo "FATAL ERROR: invalid use of -r. Please contact kgt@lanl.gov."
    print_use; exit 1
  fi
  # special defaults for regress_mode
  if [[ -d /scratch/regress/logs ]]; then
    # ccs-net machines
    logdir=/scratch/regress/logs
  elif [[ -d /usr/projects/jayenne/regress/logs ]]; then
    # HPC machines, Darwin
    logdir=/usr/projects/jayenne/regress/logs
  fi
fi
export regress_mode logdir

if ! [[ -d $logdir ]]; then
  run "mkdir -p $logdir" || die "Cannot create logdir = $logdir"
fi

#------------------------------------------------------------------------------#
# Banner
#------------------------------------------------------------------------------#

echo "Starting checkpr.sh $@ ..."
echo "   project      = $project"
echo "   pr           = $pr"
echo "   target       = $target"
echo "   scratchdir   = $scratchdir"
echo "   logdir       = $logdir"
echo "   rscriptdir   = $rscriptdir"
echo "   regress_mode = $regress_mode"

#------------------------------------------------------------------------------#
# Kick off a regression
#------------------------------------------------------------------------------#

function startCI()
{
  local project=$1
  local build_type=$2
  local extra=$3
  local pr=$4
  if ! [[ ${pr} == "develop" ]]; then
    pr=pr${pr}
  fi
  if [[ ${project} == "draco" ]] && [[ ${extra} == 'vtest' ]]; then
    # Capsaicin/Jayenne -e vtest uses Draco w/o 'vtest'
    extra="na"
  fi
  if [[ ${extra} == 'na' ]]; then
    extra=""
    edash=""
    eflag=""
  elif [[ ${extra} == 'autodoc' ]]; then
    extra=""
    edash=""
    eflag=""
    autodoc='-a'
  else
    extrastring="(${extra}) "
    edash="-"
    eflag="-e"
  fi
  if [[ ${regress_mode} == "on" ]]; then
    rflag="-r"
  fi

  # logfile=${logdir}/${machine_name_short}-${project}-${build_type}${edash}${extra}-master-${pr}.log
  # allow_file_to_age $logfile 600
  echo -e "\n- Starting CI regression ${extrastring}for ${pr}."
#  echo "  Log: $logfile"
  echo "  Log: $logdir/${machine_name_short}-${build_type}-master-YYYYMMDD-hhmm.log"
  echo " "
  cmd="$rscriptdir/regression-master.sh ${rflag} -b ${build_type}"
  cmd="$cmd ${autodoc} ${eflag} ${extra} -p ${project} -f ${pr}"
  case $target in
    ccscs* )  # build one at a time.
      ;;
    * )
      # Run all builds simultaneously (via job submission system)
      if ! [[ $project == "draco" ]]; then
        cmd="$cmd &"
      fi
      ;;
  esac
  echo "$cmd"
  eval "$cmd"
}

#------------------------------------------------------------------------------#
# Prepare for Jayenne/Capsaicin
# Do we need to build/install a new draco-develop to link against?
#------------------------------------------------------------------------------#

case $target in
  ccscs*) machine_name_short=ccscs ;;
  sn*) machine_name_short=sn ;;
  tt*) machine_name_short=tt ;;
  darwin*) machine_name_short=darwin ;;
esac

# This file tracks when draco was last built.
# reset the draco-last-built tag (draco has changed)
draco_tag_file=$logdir/last-draco-develop-${target}.log
if [[ $rmlastdracotag == "on" ]]; then
  if [[ -f $draco_tag_file ]]; then
    run "rm $draco_tag_file"
  fi
fi

case $project in
  jayenne|capsaicin)
    # Do we need to build draco? Only build draco-develop once per day.
    eval "$(date +'today=%F now=%s')"
    midnight=$(date -d "$today 0" +%s)
    if [[ -f $draco_tag_file ]]; then
      draco_last_built=$(date +%s -r $draco_tag_file)
    else
      draco_last_built=0
    fi
    if [[ $midnight -gt $draco_last_built ]]; then
      echo " "
      echo "Found a Jayenne or Capsaicin PR, but we need to build draco-develop first..."
      echo " "

      # Call this script recursively to build the draco 'develop' branch.
      if [[ ${regress_mode} == "on" ]]; then
        rflag="-r"
      fi
      run "$rscriptdir/checkpr.sh ${rflag} -p draco"
      echo " "

      # Reset the modified date on the file used to determine when draco was
      # last built.
      echo "date &> $draco_tag_file"
      date &> $draco_tag_file
    fi
    ;;
esac

#------------------------------------------------------------------------------#
# Process the build.
#------------------------------------------------------------------------------#

echo " "
case $target in

  # CCS-NET: Release, vtest, coverage
  ccscs2*)
    startCI ${project} Release autodoc $pr
    startCI ${project} Debug coverage $pr
    ;;

  # CCS-NET: Release, vtest, coverage
  ccscs3*)
    if [[ ${project} == "draco" ]]; then
      startCI ${project} Release na $pr
    else
      startCI ${project} Release vtest $pr
    fi
    startCI ${project} Debug clang $pr
    ;;

  # CCS-NET: Valgrind (Debug)
  ccscs6*) startCI ${project} Debug valgrind $pr ;;

  # Snow: Debug
  sn-fe*)
    startCI ${project} Release na $pr
    if ! [[ ${project} == "draco" ]]; then
      startCI ${project} Release vtest $pr
    fi
    startCI ${project} Debug fulldiagnostics $pr
    ;;

  # Trinitite: Release
  tt-fe*)
    startCI ${project} Release na $pr
    if ! [[ ${project} == "draco" ]]; then
      startCI ${project} Release vtest $pr;
    fi
    startCI ${project} Release knl $pr
    ;;

  # Darwin: Disabled
  darwin-fe*)
    # startCI ${project} Release na $pr
    ;;

  # These cases are not automated checks of PRs.  However, these machines are
  # supported if this script is started by a developer:
  ccscs[14]*)
    startCI ${project} Release na $pr
    startCI ${project} Debug na $pr ;;
  ccscs[589]*)
    startCI ${project} Debug coverage $pr ;;
  ba-fe* | pi-fe* | wf-fe*)
    startCI ${project} Release na $pr
    startCI ${project} Release vtest $pr
    startCI ${project} Debug fulldiagnostics $pr
    ;;

  *) echo "Unknown target machine: target = $target" ; exit 1 ;;

esac

# Wait for all subprocesses to finish before exiting this script
if [[ `jobs -p | wc -l` -gt 0 ]]; then
  echo " "
  echo "Jobs still running (if any):"
  for job in `jobs -p`; do
    echo "  waiting for job $job to finish..."
    wait $job
    echo "  waiting for job $job to finish...done"
  done
fi

echo " "
echo "All done."

#------------------------------------------------------------------------------#
# end checkpr.sh
#------------------------------------------------------------------------------#
