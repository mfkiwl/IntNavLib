#!/bin/sh

set -e

REPO_ROOT="$(git rev-parse --show-toplevel)"

cd "$REPO_ROOT/build"
ctest --output-on-failure --verbose
cd "$REPO_ROOT"