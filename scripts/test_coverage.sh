#!/bin/sh

set -e

REPO_ROOT="$(git rev-parse --show-toplevel)"

cd "$REPO_ROOT/scripts"

# build in debug mode, with coverage flag
./build.sh --coverage Debug

cd "$REPO_ROOT/build"

# launch tests -> creates coverage files
ctest

# generate coverage report
lcov --capture --directory . --output-file coverage.info
# remove system libraries from check
lcov --remove coverage.info '/usr/*' -o coverage.info
genhtml coverage.info --output-directory coverage_html
firefox coverage_html/index.html

cd "$REPO_ROOT"