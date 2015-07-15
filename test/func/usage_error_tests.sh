#!/bin/bash

. test_functions.sh

### Wrong commands
run $GQT foo
assert_fail_to_stderr $EX_USAGE $LINENO

### Missing command line attributes
run $GQT 
assert_fail_to_stderr $EX_USAGE $LINENO

for CMD in pca-shared gst fst calpha query "convert ped" "convert bcf"
do
    run $GQT $CMD 
    assert_fail_to_stderr $EX_USAGE $LINENO
done

### Missing file
for CMD in pca-shared gst fst calpha query "convert ped" "convert bcf"
do
    run $GQT $CMD -i no_such_file
    assert_fail_to_stderr $EX_NOINPUT $LINENO
done
