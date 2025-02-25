#!/bin/bash

source shared.sh

all_passed=1

function run_single_test {
    echo "Running test $1..."
    rerun_local_iokerneld
    if [[ $1 == *"tcp"* ]]; then
    	rerun_mem_server
    fi
    run_program ./bin/$1
}

function run_all_tests {
    TESTS=`ls bin | grep test_`
    for test in $TESTS
    do
        run_single_test $test
    done
}

function cleanup {
    kill_local_iokerneld
    kill_mem_server
}

# run_all_tests
run_single_test $*

cleanup

if [[ $all_passed -eq 1 ]]; then
    exit 0
else
    exit -1
fi
