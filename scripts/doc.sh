#!/bin/sh

REPO_ROOT="$(git rev-parse --show-toplevel)"

cd "$REPO_ROOT"
doxygen Doxyfile