#!/bin/bash
##---------------------------------------------------------------------------##
## File  : regression/repository_list.sh
## Date  : Thursday, May 03, 2018, 08:47 am
## Author: Kelly Thompson
## Note  : Copyright (C) 2018, Triad National Security, LLC.
##         All rights are reserved.
##---------------------------------------------------------------------------##

#------------------------------------------------------------------------------#
# A list of repositories that the following scripts will use:
# 1. sync_repositories.sh - make a 'bare' copy of each listed repository and
#    store at ccscs2:/ccs/codes/radtran/git; every 20 minutes update repo
#    looking for new PR branches; start CI if new PR branch found.
# 2. push_repositories_xf.sh - for each listed repository, tar and push to the
#    red network.
# 3. pull_repositories_xf.sh - for each listed repository, pull and expand
#    archive.  Repositories will appear on the red network at
#    /usr/projects/draco/git/.
#------------------------------------------------------------------------------#

# Github.com/lanl repositories:

# DRACO: For all machines running this script, copy all of the git repositories
# to the local file system.
github_projects=(
  lanl/Draco
  lanl/branson
)

# Gitlab.lanl.gov repositories:

# JAYENNE, CAPSAICIN, CSK: For all machines running this scirpt, copy all of the
# git repositories to the local file system.
gitlab_projects=(
  CSK/CSK
  Draco/dracodoc
  capsaicin/capsaicin
  capsaicin/core
  capsaicin/docs
  capsaicin/npt
  capsaicin/trt
  jayenne/imcdoc
  jayenne/jayenne
  user_contrib/user_contrib
)

#------------------------------------------------------------------------------#
# Generate a combined list
git_projects+=("${github_projects[@]}")
git_projects+=("${gitlab_projects[@]}")

#------------------------------------------------------------------------------#
# End repository_list.sh
#------------------------------------------------------------------------------#
