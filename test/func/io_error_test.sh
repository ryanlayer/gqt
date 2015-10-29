#!/bin/bash

. test_functions.sh

### Missing file
for CMD in pca-shared gst fst calpha query "convert ped" "convert bcf"
do
    run $GQT $CMD -i no_such_file.gqt
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

run $GQT query -i $BCF.gqt -V no_file
assert_fail_to_stderr $EX_NOINPUT $LINENO

run $GQT query -i $BCF.gqt -B no_file
assert_fail_to_stderr $EX_NOINPUT $LINENO

# corrupt GQT file
make_index
echo 11111111111111111111111 > $BCF.gqt
run $GQT query -i $BCF.gqt 
assert_fail_to_stderr $EX_IOERR $LINENO
assert_in_stderr $BCF.gqt $LINENO
rm_index

# corrupt BIM file
make_index
echo 11111111111111111111111 > $BCF.bim
run $GQT query -i $BCF.gqt 
assert_fail_to_stderr $EX_IOERR $LINENO
assert_in_stderr $BCF.bim $LINENO
rm_index

# corrupt VID file
make_index
echo 11111111111111111111111 > $BCF.vid
run $GQT query -i $BCF.gqt 
assert_fail_to_stderr $EX_IOERR $LINENO
assert_in_stderr $BCF.vid $LINENO
rm_index

# corrupt VID file
make_index
perl -e '$out = "foobar";print pack("V*",10,0,length($out),0,48,0,1,0,666,666);print $out' > truncated.bin
run $GQT query -i $BCF.gqt  -B truncated.bin
assert_fail_to_stderr $EX_IOERR $LINENO
assert_in_stderr truncated.bin $LINENO
rm_index
rm truncated.bin

# empty PED file
make_index
echo "" > empty.ped
run $GQT convert ped -i $BCF  -p empty.ped
assert_fail_to_stderr $EX_NOINPUT $LINENO
assert_in_stderr "gqt: Empty PED file 'empty.ped'." $LINENO
rm empty.ped
rm_index

# PED with bad field
make_index
echo -e "A\tf()" > bad_field.ped
run $GQT convert ped -i $BCF  -p bad_field.ped
assert_fail_to_stderr $EX_NOINPUT $LINENO
assert_in_stderr "gqt: Invalid character '(' in field name 'f()' from file 'bad_field.ped'" $LINENO
rm bad_field.ped
rm_index

# PED with fields and no data
make_index
echo -e "A\tB" > bad_field.ped
run $GQT convert ped -i $BCF  -p bad_field.ped
assert_fail_to_stderr $EX_NOINPUT $LINENO
assert_in_stderr "No data in PED file 'bad_field.ped'" $LINENO
#rm bad_field.ped
rm_index

# PED with short line
make_index
echo -e "A\tB\n1\t2\n1\n" > bad_field.ped
run $GQT convert ped -i $BCF  -p bad_field.ped
assert_fail_to_stderr $EX_NOINPUT $LINENO
assert_in_stderr "Missing field in file 'bad_field.ped' on line 3." $LINENO
rm bad_field.ped
rm_index

# PED with long line
make_index
echo -e "A\tB\n1\t2\n1\t2\t3\n\n" > bad_field.ped
run $GQT convert ped -i $BCF  -p bad_field.ped
assert_fail_to_stderr $EX_NOINPUT $LINENO
assert_in_stderr "Extra field in file 'bad_field.ped' on line 3." $LINENO
rm bad_field.ped
rm_index

# PED file extra field in one row
make_index
echo -e "A\tB\n1\t2\n1\t2" > bad_field.ped
run $GQT convert ped -i $BCF  -p bad_field.ped -c 3
assert_fail_to_stderr $EX_NOINPUT $LINENO
assert_in_stderr "gqt: Too few columns in PED file 'bad_field.ped'. Sample IDs column specified as 3,  but only 2 columns present." $LINENO
rm bad_field.ped
rm_index

# PED file with extra record
make_index
echo -e "A\tB\tC" > test.ped
echo -e "1\tI0\t1" >> test.ped
echo -e "1\tI1\tx" >> test.ped
echo -e "1\tI11\t1" >> test.ped
run $GQT convert ped -i $BCF -p test.ped
assert_in_stderr "gqt: WARNING: No match found for sample 'I11' from PED file." $LINENO
rm_index

# PED file with no matching records
make_index
echo -e "A\tB\tC" > test.ped
echo -e "1\tI20\t1" >> test.ped
echo -e "1\tI30\tx" >> test.ped
echo -e "1\tI11\t1" >> test.ped
run $GQT convert ped -i $BCF -p test.ped
assert_in_stderr "gqt: WARNING: None of the samples names from column 2 in PED file test.ped matched sample names in VCF/BCF ../data/10.1e4.var.bcf'" $LINENO
assert_in_stderr "0 of 3 PED samples matched VCF/BCF database records." $LINENO
rm_index

# Test VID without VID header
make_index
run $GQT query -i $BCF.gqt -V $BCF.gqt
assert_fail_to_stderr $EX_NOINPUT $LINENO
assert_in_stderr "gqt: File '../data/10.1e4.var.bcf.gqt' is not a VID file." $LINENO
rm_index


# Test VID with VID header
make_index
run $GQT query -i $BCF.gqt -V $BCF.vid
assert_exit_code $EX_OK $LINENO
rm_index

# Test BIM without BIM header
make_index
run $GQT query -i $BCF.gqt -B $BCF.gqt
assert_fail_to_stderr $EX_NOINPUT $LINENO
assert_in_stderr "gqt: File '../data/10.1e4.var.bcf.gqt' is not a BIM file." $LINENO
rm_index

# Test BIM with BIM header
make_index
run $GQT query -i $BCF.gqt -B $BCF.bim
assert_exit_code $EX_OK $LINENO
rm_index

# Test WAHBM without WAHBM header
make_index
run $GQT query -i $BCF -G $BCF -B $BCF.bim -V $BCF.vid -d $BCF.db
assert_fail_to_stderr $EX_NOINPUT $LINENO
assert_in_stderr "gqt: File '../data/10.1e4.var.bcf' is not a GQT file." $LINENO
rm_index

# Test WAHBM with WAHBM header
make_index
run $GQT query -i $BCF -G $BCF.gqt -B $BCF.bim -V $BCF.vid -d $BCF.db
assert_exit_code $EX_OK $LINENO
rm_index

