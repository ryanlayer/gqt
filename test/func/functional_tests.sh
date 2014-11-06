#!/bin/bash

GTQ_PATH=../../bin
DATA_PATH=../data
MORE_DATA_PATH=~/data/genotq/sim/


$GTQ_PATH/gqt convert vcf-plt -f 10 -r 43 -i $DATA_PATH/10.1e4.var.vcf -o tmp.var.plt
$GTQ_PATH/gqt convert plt-ubin -i $DATA_PATH/10.1e4.ind.txt -o tmp.ubin
$GTQ_PATH/gqt convert ubin-wahbm -i tmp.ubin -o tmp.wahbm
$GTQ_PATH/gqt convert ubin-wah -i tmp.ubin -o tmp.wah
$GTQ_PATH/gqt convert ubin-plt -i tmp.ubin -o tmp.ubin.to.plt

$GTQ_PATH/gqt view plt -i $DATA_PATH/10.1e4.ind.txt > tmp.plt.plt
$GTQ_PATH/gqt view ubin -i tmp.ubin > tmp.ubin.plt
$GTQ_PATH/gqt view wahbm -i tmp.wahbm > tmp.wahbm.plt
$GTQ_PATH/gqt view wah -i tmp.wah > tmp.wah.plt

diff tmp.var.plt $DATA_PATH/10.1e4.var.txt
diff tmp.plt.plt tmp.ubin.plt
diff tmp.plt.plt tmp.wahbm.plt
diff tmp.plt.plt tmp.wah.plt
diff tmp.ubin.to.plt $DATA_PATH/10.1e4.ind.txt

rm tmp.var.plt \
    tmp.plt.plt \
    tmp.ubin.plt \
    tmp.wahbm.plt \
    tmp.wah.plt \
    tmp.ubin \
    tmp.wahbm \
    tmp.wah \
    tmp.ubin.to.plt

ARGS="-q 0 -n 5 -r 1,2,4,5,7"
$GTQ_PATH/gqt gt plt \
    -i $DATA_PATH/10.1e4.ind.txt \
    -b $DATA_PATH/10.1e4.ind.bim \
    $ARGS \
    > tmp.gt.plt
$GTQ_PATH/gqt gt ubin \
    -i $DATA_PATH/10.1e4.ind.ubin \
    -b $DATA_PATH/10.1e4.ind.bim \
    $ARGS \
    > tmp.gt.ubin
$GTQ_PATH/gqt gt wahbm \
    -i $DATA_PATH/10.1e4.ind.wahbm \
    -b $DATA_PATH/10.1e4.ind.bim \
    $ARGS \
    > tmp.gt.wahbm
$GTQ_PATH/gqt gt ipwahbm \
    -i $DATA_PATH/10.1e4.ind.wahbm \
    -b $DATA_PATH/10.1e4.ind.bim \
    $ARGS \
    > tmp.gt.ipwahbm
$GTQ_PATH/gqt gt cipwahbm\
    -i $DATA_PATH/10.1e4.ind.wahbm \
    -b $DATA_PATH/10.1e4.ind.bim \
    $ARGS \
    > tmp.gt.cipwahbm
diff tmp.gt.plt tmp.gt.ubin
diff tmp.gt.plt tmp.gt.wahbm
diff tmp.gt.plt tmp.gt.ipwahbm
diff tmp.gt.plt tmp.gt.cipwahbm
rm tmp.gt.plt \
    tmp.gt.ubin \
    tmp.gt.wahbm \
    tmp.gt.ipwahbm \
    tmp.gt.cipwahbm

ARGS="-o gt -q 0 -n 5 -r 1,2,4,5,7"
$GTQ_PATH/gqt count plt \
    -i $DATA_PATH/10.1e4.ind.txt \
    -b $DATA_PATH/10.1e4.ind.bim \
    $ARGS > tmp.count.plt

$GTQ_PATH/gqt count ubin \
    -i $DATA_PATH/10.1e4.ind.ubin \
    -b $DATA_PATH/10.1e4.ind.bim \
    $ARGS > tmp.count.ubin

$GTQ_PATH/gqt count wahbm \
    -i $DATA_PATH/10.1e4.ind.wahbm \
    -b $DATA_PATH/10.1e4.ind.bim \
    $ARGS > tmp.count.wahbm

$GTQ_PATH/gqt count ipwahbm \
    -i $DATA_PATH/10.1e4.ind.wahbm \
    -b $DATA_PATH/10.1e4.ind.bim \
    $ARGS > tmp.count.ipwahbm

$GTQ_PATH/gqt count cipwahbm \
    -i $DATA_PATH/10.1e4.ind.wahbm \
    -b $DATA_PATH/10.1e4.ind.bim \
    $ARGS > tmp.count.cipwahbm

diff tmp.count.plt tmp.count.wahbm
#diff tmp.count.plt tmp.count.ubin
#diff tmp.count.plt tmp.count.ipwahbm
#diff tmp.count.plt tmp.count.cipwahbm

rm -f tmp.count.plt \
    tmp.count.ubin \
    tmp.count.wahbm \
    tmp.count.ipwahbm \
    tmp.count.cipwahbm

plink --file ../data/10.1e4.ind --freq >/dev/null

cat plink.frq \
    | grep -v "CHR" \
    | awk '{OFS="\t";
            if ($3=="2") 
                print $1,$2,$5*$6; 
            else 
                print $1,$2,(1-$5)*$6
            }' \
    > plink.frq.common

echo -en "\n" >> plink.frq.common

$GTQ_PATH/gqt sum ipwahbm \
    -i ../data/10.1e4.ind.wahbm \
    -b ../data/10.1e4.ind.bim \
    -n 10 \
    -r 0,1,2,3,4,5,6,7,8,9 \
    -u 2 \
    -l 1 \
    > gqt.out

cat gqt.out | cut -f 1,2,7 > gqt.out.common

$GTQ_PATH/gqt sum ipwahbm \
    -a \
    -i ../data/10.1e4.ind.wahbm \
    -b ../data/10.1e4.ind.bim \
    -n 10 \
    -r 0,1,2,3,4,5,6,7,8,9 \
    -u 2 \
    -l 1 \
    > gqt.out.a

cat gqt.out.a | cut -f 1,2,7 > gqt.out.a.common

if [ -n "`diff -w gqt.out.common plink.frq.common`" ]
then 
    echo "ERROR: gqt sum does not match plink"
else
    echo "SUCCESS: gqt sum matches plink"
    rm -f gqt.out.common plink.frq.common plink.frq
fi

if [ -n "`diff -w gqt.out gqt.out.a`" ]
then 
    echo "ERROR: gqt sum does not match gqt sum -a"
else
    echo "SUCCESS: gqt sum matches gqt sum -a"
    rm -f gqt.out gqt.out.a
fi

$GTQ_PATH/gqt convert vcf-plt \
    -i ../data/10.1e4.var.vcf \
    -o .tmp.var.plt \
    -r 43 \
    -f 10

$GTQ_PATH/gqt convert plt-invert-ubin \
    -i .tmp.var.plt \
    -o .tmp.ind.ubin

$GTQ_PATH/gqt convert ubin-plt \
    -i .tmp.ind.ubin \
    -o .tmp.ind.plt

$GTQ_PATH/gqt sort plt-field-freq \
    -i .tmp.ind.plt \
    -o .tmp.sort.ind.plt

$GTQ_PATH/gqt convert plt-ubin \
    -i .tmp.sort.ind.plt \
    -o .tmp.sort.ind.ubin

$GTQ_PATH/gqt convert ubin-wahbm \
    -i .tmp.sort.ind.ubin \
    -o .tmp.sort.ind.wahbm

$GTQ_PATH/gqt convert bcf-wahbm \
    -r 43 \
    -f 10 \
    -i ../data/10.1e4.var.bcf \
    -b .tmp.bcf.sort.ind.bim \
    -o .tmp.bcf.sort.ind.wahbm

$GTQ_PATH/gqt sum ipwahbm \
    -a \
    -i .tmp.bcf.sort.ind.wahbm \
    -b .tmp.bcf.sort.ind.bim \
    -n 10 \
    -r 0,1,2,3,4,5,6,7,8,9 \
    -u 2 \
    -l 1  \
    > .tmp.bcf.sort.ind.wahbm.out


$GTQ_PATH/gqt sum ipwahbm \
    -a \
    -i .tmp.sort.ind.wahbm \
    -b .tmp.bcf.sort.ind.bim \
    -n 10 \
    -r 0,1,2,3,4,5,6,7,8,9 \
    -u 2 \
    -l 1 \
    > .tmp.bcf.sort.ind.bim.out

if [ -n "`diff -w .tmp.bcf.sort.ind.wahbm.out .tmp.bcf.sort.ind.bim.out`" ]
then 
    echo "ERROR: gqt vcf...wahbm does not match bcf-wahbm"
else
    echo "SUCCESS: gqt vcf...wahbm matches bcf-wahbm"
    rm -f .tmp.var.plt \
        .tmp.ind.ubin \
        .tmp.ind.plt \
        .tmp.sort.ind.plt \
        .tmp.sort.ind.ubin \
        .tmp.sort.ind.wahbm \
        .tmp.bcf.sort.ind.bim \
        .tmp.bcf.sort.ind.wahbm \
        .tmp.bcf.sort.ind.wahbm.out \
        .tmp.bcf.sort.ind.bim.out
fi


$GTQ_PATH/gqt convert ped-db \
    -i ../data/10.1e4.var.ped \
    -o .tmp.10.1e4.var.db

ROWS=`sqlite3 .tmp.10.1e4.var.db "select * from ped;" | wc -l`

if [ $ROWS -eq 10 ]
then
    echo "SUCCESS: Correct number of rows in PED db"
else
    echo "ERROR: 10 rows expect in PED db. $ROWS found."
fi 

ROWS=`sqlite3 .tmp.10.1e4.var.db "select * from ped where Population='ESN';"| wc -l`

if [ $ROWS -eq 2 ]
then
    echo "SUCCESS: Correct number of rows from ESN populaiton in PED db"
else
    echo "ERROR: 2 rows expect in PED db. $ROWS found."
fi 


# count the number of homo_ref rows
$GTQ_PATH/gqt sum ipwahbm \
        -i ../data/10.1e4.ind.wahbm \
        -b ../data/10.1e4.var.bim \
        -d .tmp.10.1e4.var.db \
        -q "Population ='ESN'" \
        -l 0 \
        -u 0 \
> .tmp.homo_ref.count

#ESN are the 2nd and 10th columns in the vcf
GQT_BOTH_NUM=`cat .tmp.homo_ref.count | awk '$6==2' | wc -l`
GQT_ONE_NUM=`cat .tmp.homo_ref.count | awk '$6==1' | wc -l`

VCF_BOTH_NUM=`cat ../data/10.1e4.var.vcf  | cut -f 10- | tail -n+6 | cut -f2,10 | awk '$1=="0|0" && $2=="0|0"' | wc -l`
VCF_ONE_NUM=`cat ../data/10.1e4.var.vcf  | cut -f 10- | tail -n+6 | cut -f2,10 | awk '($1=="0|0" || $2=="0|0") && !($1=="0|0" && $2=="0|0")' | wc -l`

if [ $GQT_BOTH_NUM -eq $VCF_BOTH_NUM ]
then
    echo "SUCCESS: Number of HOMO_REF in both ESN match in VCF and GQT"
else
    echo "ERROR: Number of HOMO_REF in both ESN do not match in VCF($VCF_BOTH_NUM)  and GQT($GQT_BOTH_NUM)"
fi 

if [ $GQT_ONE_NUM -eq $VCF_ONE_NUM ]
then
    echo "SUCCESS: Number of HOMO_REF in only one ESN match in VCF and GQT"
else
    echo "ERROR: Number of HOMO_REF in only one ESN do not match in VCF($VCF_ONE_NUM)  and GQT($GQT_ONE_NUM)"
fi 


$GTQ_PATH/gqt query \
        -i ../data/10.1e4.ind.wahbm \
        -b ../data/10.1e4.var.bim \
        -d .tmp.10.1e4.var.db \
        -p "Population ='ESN'" \
        -g "count(HOMO_REF) >= 1"


$GTQ_PATH/gqt query \
        -i ../data/10.1e4.ind.wahbm \
        -b ../data/10.1e4.var.bim \
        -d .tmp.10.1e4.var.db \
        -p "Population ='ITU'" \
        -g "HOMO_REF" \
        -p "Population ='ESN'" \
        -g "count(HOMO_REF) >= 1" \
        -p "Population ='ITU'" \
        -g "count(HOMO_REF) > 1"
