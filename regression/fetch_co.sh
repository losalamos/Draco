#!/bin/bash

# fetch and checkout
# $1 is the path to git
# $2 is the host ('github' or 'gitlab')
# $3 is the pull request number (e.g.: 42)

GIT=$1
githost=$2
featurebranch=$3

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

log "This is fetch_co.sh. Args = $1 $2 $3"

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

# If the branch already exists in the regression source directory, remove it and
# fetch it again.
if test "$thisbranch" = "pr${featurebranch}"; then
  # run "${GIT} reset --hard origin/develop"
  run "${GIT} checkout develop"
  run "${GIT} branch -D pr${featurebranch}"
  # run "${GIT} pull origin $prdir/${featurebranch}/head"
fi
run "${GIT} fetch origin ${prdir}/${featurebranch}/head:pr${featurebranch}"
run "${GIT} checkout pr${featurebranch}"

# Another option is to edit .git/config:
# [core]
#         repositoryformatversion = 0
#         filemode = true
#         bare = false
#         logallrefupdates = true
# [remote "origin"]
#         url = git@github.com:lanl/Draco.git
#         fetch = +refs/heads/*:refs/remotes/origin/*
#         fetch = +refs/pull/*/head:refs/remotes/origin/pr/*
# [branch "develop"]
#         remote = origin
#         merge = refs/heads/develop
# so that 'git fetch origin' fetches all PRs and 'git checkout pr/999' will be a
# tracking branch. The problem with this approach is that all PRs are downloaded
# and this can take time.
