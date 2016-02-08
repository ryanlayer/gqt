#!/bin/bash

. test_functions.sh

### Wrong commands
run $GQT foo
assert_fail_to_stderr $EX_USAGE $LINENO

### Missing command line attributes
run $GQT 
assert_fail_to_stderr $EX_OK $LINENO

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

make_index
run $GQT query -i $BCF.gqt -p "" -g "f"
assert_fail_to_stderr $EX_USAGE $LINENO
assert_in_stderr "gqt: GENOTYPE SYNTAX ERROR: Invalid query 'f'.  Expected funciton (count, pct, maf) or gentoype (HET, HOM_REF, HOM_ALT, or UNKNOWN)." $LINENO
rm_index
 
make_index
run $GQT query -i $BCF.gqt -p "" -g "count(HET f)"
assert_fail_to_stderr $EX_USAGE $LINENO
assert_in_stderr "gqt: GENOTYPE SYNTAX ERROR: Invalid function parameter 'f'.  Expected HET, HOM_REF, HOM_ALT, or UNKNOWN." $LINENO
rm_index
 
make_index
run $GQT query -i $BCF.gqt -p "" -g "count(HET 1)"
assert_fail_to_stderr $EX_USAGE $LINENO
assert_in_stderr "gqt: GENOTYPE SYNTAX ERROR:Opperation (count,pct,maf) and condition (==,!=,<, etc.) expected prior to '1'" $LINENO
rm_index

make_index
run $GQT query -i $BCF.gqt -p "" -g "count(HET)<"
assert_fail_to_stderr $EX_USAGE $LINENO
assert_in_stderr "gqt: GENOTYPE SYNTAX ERROR: Missing condition value (after ==, <, etc.)" $LINENO
rm_index

make_index
run $GQT query -i $BCF.gqt -p "" -g "maf(HET)<1"
assert_fail_to_stderr $EX_USAGE $LINENO
assert_in_stderr "gqt: GENOTYPE SYNTAX ERROR: Opperation (maf) does not expect genotype (HOM_REF,HET,HOM_ALT,UNKNOWN) prior to ')'" $LINENO
rm_index

make_index
run $GQT query -i $BCF.gqt -p "" -g "maf()<1 f"
assert_fail_to_stderr $EX_USAGE $LINENO
assert_in_stderr "gqt: GENOTYPE SYNTAX ERROR: Invalid trailing value 'f'"
rm_index
