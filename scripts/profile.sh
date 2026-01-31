#!/bin/sh

set -e

if [ $# -lt 1 ]; then
    echo "Usage: $0 <EXE_PATH> [args...]"
    exit 1
fi

EXE_PATH="$1"
shift # shift off EXE_PATH, leaving only program args

echo "Running valgrind leak check on $EXE_PATH..."
valgrind --leak-check=full "$EXE_PATH" "$@"

echo "Running valgrind callgrind on $EXE_PATH..."
valgrind --tool=callgrind "$EXE_PATH" "$@"
# kcachegrind callgrind.out.<PID>

# Other useful valgrind tools: cachegrind, helgrind

echo "Done"
