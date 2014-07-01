#!/bin/bash

GTQ_PATH=~/src/genotq/bin/
DATA_PATH=~/src/genotq/test/data/


$GTQ_PATH/gtq convert plt-ubin -i $DATA_PATH/10.1e4.ind.txt -o tmp.ubin
$GTQ_PATH/gtq convert ubin-wahbm -i tmp.ubin -o tmp.wahbm
$GTQ_PATH/gtq convert ubin-wah -i tmp.ubin -o tmp.wah

$GTQ_PATH/gtq view plt -i $DATA_PATH/10.1e4.ind.txt > tmp.plt.plt
$GTQ_PATH/gtq view ubin -i tmp.ubin > tmp.ubin.plt
$GTQ_PATH/gtq view wahbm -i tmp.wahbm > tmp.wahbm.plt
$GTQ_PATH/gtq view wah -i tmp.wah > tmp.wah.plt

diff tmp.plt.plt tmp.ubin.plt
diff tmp.plt.plt tmp.wahbm.plt
diff tmp.plt.plt tmp.wah.plt

rm tmp.plt.plt tmp.ubin.plt tmp.wahbm.plt tmp.wah.plt
rm tmp.ubin tmp.wahbm tmp.wah

