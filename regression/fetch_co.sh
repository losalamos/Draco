#!/bin/bash

# fetch and checkout
# $1 is the path to git
# $2 is the host ('github' or 'gitlab')
# $3 is the pull request number (e.g.: 42)

GIT=$1
githost=$2
featurebranch=$3

echo "This is fetch_co.log"
echo "This is fetch_co.log" > fetch_co.log

function log()
{
  echo "$1"
  echo "$1" >> fetch_co.log
}
function run()
{
  echo "$1"
  echo "$1" >> fetch_co.log
  eval $1 >> fetch_co.log
}

if ! test -x $GIT; then
   log "FATAL ERROR: unable to run GIT=$GIT."
   exit 1
fi
if [[ ! ${githost} ]]; then
   log "FATAL ERROR: exected three arguments, example: fetch_co.sh <path to git> <github> <42>"
   log "     second argument should be 'github' or 'gitlab'"
   exit 1
fi
if [[ ! ${featurebranch} ]]; then
   log "FATAL ERROR: exected three arguments, example: fetch_co.sh <path to git> <github> <42>"
   log "     third argument should be a bare integer (the pull request number)."
   exit 1
fi

if test $githost = "github"; then
  prdir="pull"
else
  prdir="merge-requests"
fi

# Use different logic if the desired featurebranch is already checked out.
thisbranch=`git rev-parse --abbrev-ref HEAD`
log "thisbranch = $thisbranch"

if test "$thisbranch" = "pr${featurebranch}"; then
  run "${GIT} pull origin $prdir/${featurebranch}/head"
else
  run "${GIT} fetch origin ${prdir}/${featurebranch}/head:pr${featurebranch}"
  run "${GIT} checkout pr${featurebranch}"
fi
