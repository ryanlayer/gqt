PLINK_SRC=~/src/plink2
GQT_PATH=~/src/gqt/bin
TARGET_BCF=../data/10.1e4.var.bcf
TARGET_PED=../data/10.1e4.var.ped

$PLINK_SRC/plink \
    --make-bed \
    --bcf $TARGET_BCF \
    --out $TARGET_BCF.plink \
    --allow-extra-chr \
     2> /dev/null 1> /dev/null

~/src/plink2/plink \
     --bfile $TARGET_BCF.plink \
     --freq \
     --allow-extra-chr \
     --out plink.out \
     2> /dev/null 1> /dev/null

PLINK_COUNT=`cat plink.out.frq \
    | grep -v "CHR" \
    | awk '{if ($4 == "A") print $5*$6; else print (1-$5)*$6;}' \
    | awk '{s+=$1} END {print s}'`

export BCFTOOLS_PLUGINS="/Users/rl6sf/src/bcftools/plugins/"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/Users/rl6sf/src/htslib"
export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:/Users/rl6sf/src/htslib"

BCFTOOLS_COUNT=`bcftools view $TARGET_BCF\
    | bcftools plugin fill-AN-AC \
    | grep -v "#" \
    | cut -f 8 | cut -d ";" -f3 | cut -d "=" -f2 \
    | awk '{s+=$1} END {print s}'`

HET_COUNT=`$GQT_PATH/gqt query \
    -i $TARGET_BCF.gqt \
    -d $TARGET_PED.db \
    -p "" \
    -g "count(HET)" \
    | grep -v "#" \
    | cut -d "=" -f3 \
    | awk '{s+=$1} END {print s}'`


ALT_COUNT=`$GQT_PATH/gqt query \
    -i $TARGET_BCF.gqt \
    -d $TARGET_PED.db \
    -p "" \
    -g "count(HOMO_ALT)" \
    | grep -v "#" \
    | cut -d "=" -f3 \
    | awk '{s+=$1*2} END {print s}'`

GQT_COUNT=`echo $HET_COUNT $ALT_COUNT | awk '{print $1+$2;}'`

if [ $GQT_COUNT -eq $BCFTOOLS_COUNT ]
then
    echo "SUCCESS: GQT count matches BCFTOOLS count"
else
    echo "FAILURE: GQT count does not matche BCFTOOLS count. $GQT_COUNT vs $BCFTOOLS_COUNT"
fi

if [ $GQT_COUNT -eq $PLINK_COUNT ]
then
    echo "SUCCESS: GQT count matches PLINK count"
else
    echo "FAILURE: GQT count does not matche PLINK count. $GQT_COUNT vs $PLINK_COUNT"
fi

rm plink.out.frq    plink.out.log   plink.out.nosex
