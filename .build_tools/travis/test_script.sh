#!/bin/bash
# This script is meant to be called by the "script" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD

set -e

get_data() {
    # Hopefully this gets cached
    pushd data
    make
    popd
}


run_tests() {
    # first run the actual tests
    pushd bin/moanin
    make test
    popd
    pushd scripts
    make all
    popd
}

# Start by installing
pushd bin/moanin
make install
popd

# Now get data and run the tests and build the manuscript
get_data
run_tests
mkdir -p scripts/reports
