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


ARGS="-q 0 -n 5 -r 1,2,4,5,7"
$GTQ_PATH/gtq gt plt \
    -i $DATA_PATH/10.1e4.ind.txt \
    $ARGS \
    > tmp.gt.plt
$GTQ_PATH/gtq gt ubin \
    -i $DATA_PATH/10.1e4.ind.ubin \
    $ARGS \
    > tmp.gt.ubin
$GTQ_PATH/gtq gt wahbm \
    -i $DATA_PATH/10.1e4.ind.wahbm \
    $ARGS \
    > tmp.gt.wahbm
$GTQ_PATH/gtq gt ipwahbm \
    -i $DATA_PATH/10.1e4.ind.wahbm \
    $ARGS \
    > tmp.gt.ipwahbm
$GTQ_PATH/gtq gt cipwahbm\
    -i $DATA_PATH/10.1e4.ind.wahbm \
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
$GTQ_PATH/gtq count plt \
    -i $DATA_PATH/10.1e4.ind.txt \
    $ARGS > tmp.count.plt

$GTQ_PATH/gtq count ubin \
    -i $DATA_PATH/10.1e4.ind.ubin \
    $ARGS > tmp.count.ubin

$GTQ_PATH/gtq count wahbm \
    -i $DATA_PATH/10.1e4.ind.wahbm \
    $ARGS > tmp.count.wahbm

$GTQ_PATH/gtq count ipwahbm \
    -i $DATA_PATH/10.1e4.ind.wahbm \
    $ARGS > tmp.count.ipwahbm

$GTQ_PATH/gtq count cipwahbm \
    -i $DATA_PATH/10.1e4.ind.wahbm \
    $ARGS > tmp.count.cipwahbm

diff tmp.count.plt tmp.count.wahbm
diff tmp.count.plt tmp.count.ubin
diff tmp.count.plt tmp.count.ipwahbm
diff tmp.count.plt tmp.count.cipwahbm

rm tmp.count.plt \
    tmp.count.ubin \
    tmp.count.wahbm \
    tmp.count.ipwahbm \
    tmp.count.cipwahbm


$GTQ_PATH/gtq convert plt-invert \
    -i $DATA_PATH/10.1e4.ind.txt \
    -o tmp.invert 

diff $DATA_PATH/10.1e4.var.txt \
    tmp.invert 

rm tmp.invert

$GTQ_PATH/gtq convert plt-vcf \
    -i $DATA_PATH/10.1e4.var.txt \
    -o tmp.vcf 

diff tmp.vcf $DATA_PATH/10.1e4.var.vcf

rm tmp.vcf


$GTQ_PATH/gtq convert  plt-invert-ubin \
    -i $DATA_PATH/10.1e4.ind.txt \
    -o tmp.i.ubin

$GTQ_PATH/gtq view ubin -i tmp.i.ubin \
    > tmp.o.i.ubin

$GTQ_PATH/gtq view ubin -i $DATA_PATH/10.1e4.var.ubin \
    > tmp.o.var.ubin

diff tmp.o.var.ubin tmp.o.i.ubin
rm tmp.i.ubin tmp.o.i.ubin tmp.o.var.ubin
