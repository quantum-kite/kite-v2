#!/bin/bash

# This script will compare the .h5 files of all the tests
printall=0

# Each test runs from its own directory and reaches the binaries through
# $KITEX / $KITE_TOOLS, falling back to ../../build/ when run on its own.
# Resolving them here, once and absolutely, means a test still finds them if it
# is invoked from somewhere else, and that a missing build says so plainly
# instead of every test failing on "no such file".
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
: "${KITEX:=$here/../build/KITEx}"
: "${KITE_TOOLS:=$here/../build/KITE-tools}"

for binary in "$KITEX" "$KITE_TOOLS"; do
    if [[ ! -x "$binary" ]]; then
        echo "Cannot find $binary -- build the project first, or set KITEX and"
        echo "KITE_TOOLS to where your build put them."
        exit 1
    fi
done
export KITEX KITE_TOOLS

echo -e "Testing \e[1mKITEx\e[0m functionalities. \e[1m\e[32mOK\e[0m means the test passed perfectly."
for i in test0*; do
    cd $i
    echo -n "  "
    cat description | tr '\n' ' '
    bash ./test.sh quick
    cd ..
done
echo -e "\n"
