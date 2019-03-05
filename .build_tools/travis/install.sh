#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD

set -e
# Download data. Should be in parallel of the installation
pushd data
make
popd

Rscript scripts/install.R

pushd bin/moanin
make install
popd
