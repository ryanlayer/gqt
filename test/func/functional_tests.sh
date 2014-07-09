#!/bin/bash

GTQ_PATH=~/src/genotq/bin/
DATA_PATH=~/src/genotq/test/data/
MORE_DATA_PATH=~/data/genotq/sim/


$GTQ_PATH/gtq convert vcf-plt -f 10 -r 43 -i $DATA_PATH/10.1e4.var.vcf -o tmp.var.plt
$GTQ_PATH/gtq convert plt-ubin -i $DATA_PATH/10.1e4.ind.txt -o tmp.ubin
$GTQ_PATH/gtq convert ubin-wahbm -i tmp.ubin -o tmp.wahbm
$GTQ_PATH/gtq convert ubin-wah -i tmp.ubin -o tmp.wah

$GTQ_PATH/gtq view plt -i $DATA_PATH/10.1e4.ind.txt > tmp.plt.plt
$GTQ_PATH/gtq view ubin -i tmp.ubin > tmp.ubin.plt
$GTQ_PATH/gtq view wahbm -i tmp.wahbm > tmp.wahbm.plt
$GTQ_PATH/gtq view wah -i tmp.wah > tmp.wah.plt

diff tmp.var.plt $DATA_PATH/10.1e4.var.txt
diff tmp.plt.plt tmp.ubin.plt
diff tmp.plt.plt tmp.wahbm.plt
diff tmp.plt.plt tmp.wah.plt

rm tmp.var.plt \
    tmp.plt.plt \
    tmp.ubin.plt \
    tmp.wahbm.plt \
    tmp.wah.plt \
    tmp.ubin \
    tmp.wahbm \
    tmp.wah


#$GTQ_PATH/gtq view plt -i $DATA_PATH/10.1e4.ind.txt  | head -n 3 | tail -n 2
#echo
#$GTQ_PATH/gtq gt plt -i $DATA_PATH/10.1e4.ind.txt -q 0 -n 5 -r 1,2,4,5,7
#$GTQ_PATH/gtq gt ubin -i $DATA_PATH/10.1e4.ind.ubin -q 0 -n 5 -r 1,2,4,5,7 
#$GTQ_PATH/gtq gt wahbm -i $DATA_PATH/10.1e4.ind.wahbm -q 0 -n 5 -r 1,2,4,5,7 

ARGS="-q 0 -n 5 -r 1,2,4,5,7 -Q"
echo -ne "plt\t"
$GTQ_PATH/gtq gt plt -i $DATA_PATH/10.1e4.ind.txt $ARGS

echo -ne "ubin\t"
$GTQ_PATH/gtq gt ubin -i $DATA_PATH/10.1e4.ind.ubin $ARGS

echo -ne "wahbm\t"
$GTQ_PATH/gtq gt wahbm -i $DATA_PATH/10.1e4.ind.wahbm $ARGS

echo -ne "ipwahbm\t"
$GTQ_PATH/gtq gt ipwahbm -i $DATA_PATH/10.1e4.ind.wahbm $ARGS


echo

SIZE=5000
echo -ne "plt\t"
$GTQ_PATH/gtq gt plt -i $MORE_DATA_PATH/$SIZE.1e8.ind.txt $ARGS

echo -ne "ubin\t"
$GTQ_PATH/gtq gt ubin -i $MORE_DATA_PATH/$SIZE.1e8.ind.ubin $ARGS

echo -ne "wahbm\t"
$GTQ_PATH/gtq gt wahbm -i $MORE_DATA_PATH/$SIZE.1e8.ind.wah $ARGS

echo -ne "ipwahbm\t"
$GTQ_PATH/gtq gt ipwahbm -i $MORE_DATA_PATH/$SIZE.1e8.ind.wah $ARGS

ARGS="-q 0 -n 5 -r 1,2,4,5,7"
$GTQ_PATH/gtq gt ipwahbm -i $MORE_DATA_PATH/$SIZE.1e8.ind.wah $ARGS > ipwahbm.o
$GTQ_PATH/gtq gt plt -i $MORE_DATA_PATH/$SIZE.1e8.ind.txt $ARGS > plt.o
diff ipwahbm.o plt.o
