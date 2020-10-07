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

install_github_dependencies() {
    # install_github keeps failing because of the Github API rate limit. So
    # we're just going to do this by hand…
    yes | apt install pandoc-citeproc
    pushd /tmp
    wget -P . https://github.com/NelleV/timecoursedata/archive/master.zip
    unzip /tmp/master.zip
    cd timecoursedata-master
    make install-extra
    make install
    popd
    rm -rf /tmp/master.zip

    pushd /tmp
    wget -P . https://github.com/NelleV/moanin/archive/master.zip
    unzip /tmp/master.zip
    cd moanin-master
    make install-extra
    make install
    popd

    Rscript scripts/install.R
}

run_tests() {
    # first run the actual tests
    pushd scripts
    make all
    popd
}


# Now get data and run the tests and build the manuscript
# get_data
install_github_dependencies
run_tests
mkdir -p scripts/reports
