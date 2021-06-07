#!/bin/bash -l
#--------------------------------------------------------------------------------------------------#
# File  : ./.gitlab/ci/gitlab-style-checks.sh
# Date  : Tuesday, Jun 02, 2020, 12:28 pm
# Author: Kelly Thompson <kgt@lanl.gov>
# Note  : Copyright (C) 2020, Triad National Security, LLC., All rights are reserved.
#--------------------------------------------------------------------------------------------------#

# preliminaries and environment
set -e
#shellcheck source=.gitlab/ci/common.sh
source "${JAYENNE_SOURCE_DIR}/.gitlab/ci/common.sh"
#shellcheck source=.gitlab/ci/environments.sh
source "${JAYENNE_SOURCE_DIR}/.gitlab/ci/environments.sh"

echo -e "\n========== printenv ==========\n"
if ! [[ "${SLURM_NODELINST:-notset}" == "notset" ]]; then
  echo "SLURM_NODELIST = ${SLURM_NODELIST}"
fi
echo "HOSTNAME       = ${HOSTNAME}"
# printenv
# run "module avail"
printenv | grep "CI_"
echo " "

#--------------------------------------------------------------------------------------------------#
# Style check
#--------------------------------------------------------------------------------------------------#

run "cd $JAYENNE_SOURCE_DIR"
if [[ $(git branch | grep -c develop) == 0 ]]; then
  run "git branch develop origin/develop"
fi
run "${JAYENNE_SOURCE_DIR}/.gitlab/ci/check_style.sh" -t

echo -e "\n======== end .gitlab-ci-run-tests.sh ==========\n"

#--------------------------------------------------------------------------------------------------#
# End .gitlab/ci/gitlab-ci-run-tests.sh
#--------------------------------------------------------------------------------------------------#
