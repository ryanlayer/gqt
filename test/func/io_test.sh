#!/bin/bash

. test_functions.sh

# Test to see if verison and command are in the VCF header
make_index
run $GQT query -i $BCF.gqt -p "BCF_Sample in ('I0', 'I1')" -g "count(HET)"
assert_in_stdout "##gqt_queryVersion=" $LINENO
assert_in_stdout \
    "##gqt_queryCommand=$GQT query -i $BCF.gqt -p \"BCF_Sample in ('I0', 'I1')\" -g \"count(HET)\"" \
    $LINENO

run $GQT calpha -i $BCF.gqt \
            -p "BCF_Sample in ( 'I0', 'I1', 'I2', 'I3', 'I4')"\
            -p "BCF_Sample in ( 'I5', 'I6', 'I7', 'I8', 'I9')" 
assert_in_stdout "##gqt_calphaVersion=" $LINENO
assert_in_stdout \
    "##gqt_calphaCommand=$GQT calpha -i $BCF.gqt -p \"BCF_Sample in ( 'I0', 'I1', 'I2', 'I3', 'I4')\" -p \"BCF_Sample in ( 'I5', 'I6', 'I7', 'I8', 'I9')\"" \
    $LINENO 

run $GQT fst -i $BCF.gqt -p "BCF_Sample in ( 'I0', 'I1', 'I2', 'I3', 'I4')"\
    -p "BCF_Sample in ( 'I5', 'I6', 'I7', 'I8', 'I9')" 
assert_in_stdout "##gqt_fstVersion=" $LINENO
assert_in_stdout \
    "##gqt_fstCommand=$GQT fst -i $BCF.gqt -p \"BCF_Sample in ( 'I0', 'I1', 'I2', 'I3', 'I4')\" -p \"BCF_Sample in ( 'I5', 'I6', 'I7', 'I8', 'I9')\"" \
    $LINENO


run $GQT gst -i $BCF.gqt -p "BCF_Sample in ( 'I0', 'I1', 'I2', 'I3', 'I4')"\
    -p "BCF_Sample in ( 'I5', 'I6', 'I7', 'I8', 'I9')" 
assert_in_stdout "##gqt_gstVersion=" $LINENO
assert_in_stdout \
    "##gqt_gstCommand=$GQT gst -i $BCF.gqt -p \"BCF_Sample in ( 'I0', 'I1', 'I2', 'I3', 'I4')\" -p \"BCF_Sample in ( 'I5', 'I6', 'I7', 'I8', 'I9')\"" \
    $LINENO
 

make_index
echo -e "A\tB\tC" > test.ped
echo -e "1\t1\t1" >> test.ped
echo -e "1\t1\tx" >> test.ped
echo -e "1\ty\t1" >> test.ped
run $GQT convert ped -i $BCF -p test.ped
run sqlite3 test.ped.db .schema
assert_in_stdout "BCF_ID INTEGER, BCF_Sample TEXT, A INTEGER, B TEXT, C TEXT" $LINENO
rm test.ped
rm_index

make_index
echo -e "A\tB\tC" > test.ped
echo -e "1\tI0\t1" >> test.ped
echo -e "1\tI1\tx" >> test.ped
echo -e "1\tI2\t1" >> test.ped
run $GQT convert ped -i $BCF -p test.ped
assert_in_stderr "3 of 3 PED samples matched VCF/BCF database records."
rm_index
