#!/bin/bash

. test_functions.sh

### Missing file
for CMD in pca-shared gst fst calpha query "convert ped" "convert bcf"
do
    run $GQT $CMD -i no_such_file
    assert_fail_to_stderr $EX_NOINPUT $LINENO
done

### Bad BCF file
echo XYZ > bad.bcf
run $GQT convert bcf -i bad.bcf
assert_fail_to_stderr $EX_DATAERR $LINENO

# tests for bad bcf in init_bcf_file(char *file_name)
run $GQT convert bcf -i bad.bcf -r 43 -f 10
assert_fail_to_stderr $EX_DATAERR $LINENO

rm bad.bcf

rm -f $BCF.tbi $BCF.csi
run $GQT convert bcf -i $BCF
assert_fail_to_stderr $EX_NOINPUT $LINENO

$BCFTOOLS view -Ov $BCF > $BCF.gz

run $GQT convert bcf -i $BCF.gz
assert_fail_to_stderr $EX_NOINPUT $LINENO

rm $BCF.gz 

$BCFTOOLS index $BCF

run $GQT convert bcf -i $BCF -t /not_a_dir/
assert_fail_to_stderr $EX_CANTCREAT $LINENO

run $GQT convert ped -i $BCF -p no_file
assert_fail_to_stderr $EX_NOINPUT $LINENO

run $GQT convert ped -i no_file -p no_file
assert_fail_to_stderr $EX_NOINPUT $LINENO

$GQT convert bcf -i $BCF 2> /dev/null
$GQT convert ped -i $BCF 2> /dev/null

run $GQT query -i $BCF.gqt -p "Bar==1" -g "HET"
assert_fail_to_stderr $EX_SOFTWARE $LINENO

run $GQT query -i $BCF.gqt -d no_db -p "Bar==1" -g "HET"
assert_fail_to_stderr $EX_NOINPUT $LINENO

run $GQT query -i $BCF.gqt -v no_file
assert_fail_to_stderr $EX_NOINPUT $LINENO

run $GQT query -i $BCF.gqt -b no_file
assert_fail_to_stderr $EX_NOINPUT $LINENO


