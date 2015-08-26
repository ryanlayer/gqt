#!/bin/bash

BCFTOOLS=bcftools
VCFTOOLS=vcftools
PLINK=plink
GQT=../../bin/gqt
SQLITE=sqlite3
DATA_PATH=../data

if [[ ! -f "`which $BCFTOOLS`" ]]
then
    echo -e "ERROR($LINENO): BCFTOOLS not found.  Please set path in functional_tests.sh"
    exit 1
elif [[ ! -f "`which $PLINK`" ]]
then
    echo -e "ERROR($PLINK): PLINK not found.  Please set path in functional_tests.sh"
elif [[ ! -f "`which $GQT`" ]]
then
    echo -e "ERROR($GQT): PLINK not found.  Please set path in functional_tests.sh"
    exit 1
elif [[ ! -f "`which $SQLITE`" ]]
then
    echo -e "ERROR($SQLITE): SQLITE3 not found.  Please set path in functional_tests.sh"
fi


BCF=$DATA_PATH/10.1e4.var.bcf

clean_up()
{
    ls $DATA_PATH/* \
        | grep -v "^$BCF$" \
        | grep -v "^$DATA_PATH/diff_gts.bcf$" \
        | grep -v "^$DATA_PATH/more_fields.ped$" \
        | grep -v "^$DATA_PATH/too_many_fields.ped$" \
        | xargs rm
}

clean_up

$GQT convert bcf \
    -i $BCF \
    2> /dev/null

if [ $? -ne 0 ]
then
    echo "SUCCESS($LINENO): Fail to autodetect number records fields without index"
else
    echo "ERROR($LINENO): Did not fail to autodetect number records fields without index"
    exit
fi

$BCFTOOLS index $BCF

$GQT convert bcf \
    -i $BCF \
    2> /dev/null

if [ $? -eq 0 ]
then
    echo "SUCCESS($LINENO): Autodetect number records fields with index"
else
    echo "ERROR($LINENO): Fail to autodetect number records fields with index"
    exit
fi

if [ -e "$BCF.gqt" ]
then
    echo "SUCCESS($LINENO): Autocomplete GQT index name from BCF nam"
else
    echo "ERROR($LINENO): Did not autocomplete GQT index name from BCF name"
    exit
fi

if [ -e "$BCF.vid" ]
then
    echo "SUCCESS($LINENO): Autocomplete VID index name from BCF nam"
else
    echo "ERROR($LINENO): Did not autocomplete VID index name from BCF name"
    exit
fi

if [ -e "$BCF.bim" ]
then
    echo "SUCCESS($LINENO): Autocomplete BIM index name from BCF nam"
else
    echo "ERROR($LINENO): Did not autocomplete BIM index name from BCF name"
    exit
fi

rm -f tmp.bcf.bim \
      tmp.bcf.vid \
      tmp.bcf.gqt 
 
$GQT convert bcf \
    -r 43 \
    -f 10 \
    -i $BCF \
    -b tmp.bcf.bim \
    -v tmp.bcf.vid \
    -o tmp.bcf.gqt \
    2> /dev/null

$GQT convert ped \
    -i $BCF \
    -o tmp.bcf.db \
    2> /dev/null

$GQT convert bcf \
    -r 43 \
    -f 10 \
    -i $BCF \
    2> /dev/null

$GQT convert ped \
    -i $BCF \
    2> /dev/null

$GQT query -i tmp.bcf.gqt -b tmp.bcf.bim -v tmp.bcf.vid -d tmp.bcf.db\
| grep -v "gqt_queryCommand" > tmp.spec

$GQT query -i $BCF.gqt \
| grep -v "gqt_queryCommand" > tmp.nospec

if diff tmp.spec tmp.nospec > /dev/null
then
    echo "SUCCESS($LINENO): Auto output BIM, VID, PED  matches specified"
    rm tmp.spec tmp.nospec
else
    echo "ERROR($LINENO): Auto output BIM, VID, PED does not matche specified "
    exit
fi

$GQT convert ped \
    -i $BCF 2>/dev/null

if [ -e "$BCF.db" ]
then
    echo "SUCCESS($LINENO): Auto output file on ped convert correct"
else
    echo "ERRROR: Auto output file on ped convert not correct"
fi

rm -f tmp.ped.db

$GQT convert ped \
    -i $BCF \
    -o tmp.bcf.db \
    2> /dev/null

if [ -e "tmp.bcf.db" ]
then
    echo "SUCCESS($LINENO): Specified output file on ped convert correct"
else
    echo "ERRROR: Specified output file on ped convert not correct"
fi

if diff tmp.bcf.db $BCF.db > /dev/null
then
    echo "SUCCESS($LINENO): Auto output PED DB matches specified PED DB"
    rm tmp.bcf.db
else
    echo "ERROR($LINENO): Auto output PED DB does not match specified BED DB"
    exit
fi

if [[ -f "`which $SQLITE`" ]]
then
    ROWS=`$SQLITE $BCF.db "select * from ped;" | wc -l`

    if [ $ROWS -eq 10 ]
    then
        echo "SUCCESS($LINENO): Correct number of rows in PED db"
    else
        echo "ERROR($LINENO): 10 rows expect in PED db. $ROWS found."
        exit
    fi 
else
    echo "SKIP($LINENO): SQLITE3 not set"
fi

$GQT convert ped \
    -i $BCF \
    -p $DATA_PATH/more_fields.ped \
    2> /dev/null

if [ -e "$DATA_PATH/more_fields.ped.db" ]
then
    echo "SUCCESS($LINENO): Specified output file on ped convert correct"
else
    echo "ERRROR: Specified output file on ped convert not correct"
fi

if [[ -f "`which $SQLITE`" ]]
then
    ROWS=`$SQLITE $DATA_PATH/more_fields.ped.db "select * from ped WHERE BCF_Sample == Individual_ID;" | wc -l`

    if [ $ROWS -eq 10 ]
    then
        echo "SUCCESS($LINENO): Correct number of rows in PED db"
    else
        echo "ERROR($LINENO): 10 rows expect in PED db. $ROWS found."
        $SQLITE $DATA_PATH/more_fields.ped.db "select * from ped"
        exit
    fi 
else
    echo "SKIP($LINENO): SQLITE3 not set"
fi

$GQT convert ped \
    -i $BCF \
    -p $DATA_PATH/too_many_fields.ped \
    2> /dev/null

if [[ -f "`which $SQLITE`" ]]
then
    ROWS=`$SQLITE $DATA_PATH/too_many_fields.ped.db "select * from ped;" | wc -l`
    if [ $ROWS -eq 10 ]
    then
        echo "SUCCESS($LINENO): Correct number of rows in PED db when PED has extra rows"
    else
        echo "ERROR($LINENO): 10 rows expect in PED db. $ROWS found."
        exit
    fi 
else
    echo "SKIP($LINENO): SQLITE3 not set"
fi


# count the number of homo_ref rows
GQT_BOTH_NUM=`$GQT query \
    -i $BCF.gqt \
    -d $DATA_PATH/more_fields.ped.db \
    -p "Population ='ESN'" \
    -g "HOM_REF" \
    | grep -v "#" \
    | wc -l`

S=`cat $DATA_PATH/more_fields.ped \
    | awk '$7=="ESN"' \
    | cut -f 2 \
    | paste -d, - -`

VCF_BOTH_NUM=`$BCFTOOLS view -s "$S" $BCF \
    | grep -v "^#" \
    | awk '$10 =="0|0" && $11=="0|0"' \
    | wc -l`

if [ "$GQT_BOTH_NUM" -eq "$VCF_BOTH_NUM" ]
then
    echo "SUCCESS($LINENO): Number of HOM_REF in both ESN match in VCF and GQT"
else
    echo "ERROR($LINENO): Number of HOM_REF in both ESN do not match in VCF($VCF_BOTH_NUM)  and GQT($GQT_BOTH_NUM)"
    echo -e"
    $GQT query \
        -i $BCF.gqt \ 
        -d $DATA_PATH/more_fields.ped.db \ 
        -p "Population ='ESN'" \
        -g "HOM_REF" \
        | grep -v "#" \
        | wc -l"
        exit
fi 


$GQT query \
    -i $BCF.gqt \
    -d $DATA_PATH/more_fields.ped.db \
    -p "" \
    -g "pct(HOM_REF)" \
    | grep -v "^#" \
    | awk '{
        split($8,a,";");
        for (i=1;i<=length(a);i++) {
            split(a[i],b,"=");
            if (b[1]=="GQT_0")
                print b[2];
        }
    }' > tmp.pct

$GQT query \
    -i $BCF.gqt \
    -d $DATA_PATH/more_fields.ped.db \
    -p "" \
    -g "count(HOM_REF)" \
    | grep -v "^#" \
    | awk '{
        split($8,a,";");
        for (i=1;i<=length(a);i++) {
            split(a[i],b,"=");
            if (b[1]=="GQT_0")
                print b[2];
        }
    }' > tmp.count

if [ -n "`paste tmp.pct tmp.count | awk '$1*10 != $2'`" ]
then 
    L=$LINENO
    mv tmp.pct tmp.$L.pct
    mv tmp.count tmp.$L.count
    echo "ERROR($L): gqt query pct does not match gqt query count"
    exit
else
    echo "SUCCESS($LINENO): gqt query pct matches gqt query count"
    rm tmp.pct tmp.count 
fi

$GQT query \
    -i $BCF.gqt \
    -d $DATA_PATH/more_fields.ped.db \
    -p "" \
    -g "pct(HOM_REF)" \
    -p "" \
    -g "count(HOM_REF)" \
    > tmp.gqt

cat tmp.gqt \
    | grep -v "#" \
    | awk '{
        split($8,a,";");
        for (i=1;i<=length(a);i++) {
            split(a[i],b,"=");
            if (b[1]=="GQT_0")
                print b[2];
        }
    }' > tmp.pct

cat tmp.gqt \
    | grep -v "#" \
    | awk '{
        split($8,a,";");
        for (i=1;i<=length(a);i++) {
            split(a[i],b,"=");
            if (b[1]=="GQT_1")
                print b[2];
        }
    }' > tmp.count

if [ -n "`paste tmp.pct tmp.count | awk '$1*10 != $2'`" ]
then 
    L=$LINENO
    mv tmp.pct tmp.$L.pct
    mv tmp.count tmp.$L.count
    mv tmp.gqt tmp.$L.gqt
    echo "ERROR($L): gqt multi query pct and count do not match"
    exit
else
    echo "SUCCESS($LINENO): gqt multi query pct and count match" 
    rm tmp.gqt tmp.count tmp.pct
fi

$BCFTOOLS view -c 10 $BCF \
    | grep -v "^#" \
    | cut -f1-3 \
    > tmp.bcf

$GQT query \
    -i $BCF.gqt \
    -d $DATA_PATH/more_fields.ped.db \
    -p "" \
    -g "maf()>=0.5" \
    | grep -v "^#" \
    | cut -f1-3 \
    > tmp.gqt

if diff tmp.bcf tmp.gqt > /dev/null
 then 
    echo "SUCCESS($LINENO): BCFTOOLS nref count does not match GQT maf"
    rm tmp.gqt tmp.bcf
else
    L=$LINENO
    mv tmp.bcf tmp.$L.bcf
    mv tmp.gqt tmp.$L.gqt
    echo "ERROR($L): BCFTOOLS nref count does not match GQT maf"
    exit
fi


$BCFTOOLS view -Ov $BCF -o $DATA_PATH/10.1e4.var.vcf

$GQT convert bcf \
    -r 43 \
    -f 10 \
    -i $DATA_PATH/10.1e4.var.vcf \
    2> /dev/null

$BCFTOOLS view -Oz $BCF -o $DATA_PATH/10.1e4.var.vcf.gz

$GQT convert bcf \
    -r 43 \
    -f 10 \
    -i $DATA_PATH/10.1e4.var.vcf.gz\
    2> /dev/null

$GQT query \
    -i $BCF.gqt \
    -d $DATA_PATH/more_fields.ped.db \
    -p "" \
    -g "maf()>=0.5" \
    | grep -v "^#" \
    | cut -f1-3 \
    > tmp.bcf.gqt

$GQT query \
    -i $DATA_PATH/10.1e4.var.vcf.gqt \
    -d $DATA_PATH/more_fields.ped.db \
    -p "" \
    -g "maf()>=0.5" \
    | grep -v "^#" \
    | cut -f1-3 \
    > tmp.vcf.gqt

if diff tmp.bcf.gqt tmp.vcf.gqt > /dev/null
 then 
    echo "SUCCESS($LINENO): BCF-based GQT matches VCF-based GQT"
    rm tmp.vcf.gqt
else
    L=$LINENO
    cp tmp.bcf.gqt tmp.$L.bcf.gqt
    mv tmp.vcf.gqt tmp.$L.vcf.gqt
    echo "ERROR($L): BCF-based GQT does not matches VCF-based GQT"
    exit
fi

$GQT query \
    -i $DATA_PATH/10.1e4.var.vcf.gz.gqt \
    -d $DATA_PATH/more_fields.ped.db \
    -p "" \
    -g "maf()>=0.5" \
    | grep -v "^#" \
    | cut -f1-3 \
    > tmp.vcf.gz.gqt

if diff tmp.bcf.gqt tmp.vcf.gz.gqt > /dev/null
 then 
    echo "SUCCESS($LINENO): BCF-based GQT matches VCF.GZ-based GQT"
    rm tmp.vcf.gz.gqt
else
    L=$LINENO
    cp tmp.bcf.gqt tmp.$L.bcf.gqt
    mv tmp.vcf.gz.gqt tmp.$L.vcf.gz.gqt
    echo "ERROR($L): BCF-based GQT does not matches VCF.GZ.-based GQT"
    exit
fi

rm tmp.bcf.gqt


#rm -f tmp.ped
#echo -ne "Family ID\tIndividual ID\tPaternal ID\tMaternal ID\tGender\tPhenotyp\n" > tmp.ped
#echo -ne "Y025\tNA18907\t0\t0\t2\t0\n" >> tmp.ped
#echo -ne "NG108\tHG03519\tHG03518\tHG03517\t1\t0\n" >> tmp.ped
#echo -ne "m027\tNA19758\t0\t0\t2\t0\n" >> tmp.ped
#echo -ne "IT060\tHG04015\t0\t0\t1\t0\n" >> tmp.ped
#echo -ne "test space here\tand here\there too\tone more\t1\t0\n" >> tmp.ped
#
#if [[ -f "`which $SQLITE`" ]]
#then
#    $GQT convert ped \
#    -i tmp.ped
#
#    SPACE_R0=`$SQLITE tmp.ped.db "select Ind_ID from ped where Family_ID = 'test space here';"`
#    SPACE_R1=`$SQLITE tmp.ped.db "select Ind_ID from ped where Individual_ID = 'and here';"`
#    SPACE_R2=`$SQLITE tmp.ped.db "select Ind_ID from ped where Paternal_ID = 'here too';"`
#    SPACE_R3=`$SQLITE tmp.ped.db "select Ind_ID from ped where Maternal_ID = 'one more';"`
#
#    if [ $SPACE_R0 -eq $SPACE_R1 ]
#    then
#    if [ $SPACE_R0 -eq $SPACE_R2 ]
#    then
#        if [ $SPACE_R0 -eq $SPACE_R3 ]
#        then
#            echo "SUCCESS($LINENO): Spaces acceped in cell values"
#            rm tmp.ped tmp.ped.db
#        else
#            echo "ERROR($LINENO): Spaces not acceped in cell values"
#        fi
#    else
#        echo "ERROR($LINENO): Spaces not acceped in cell values"
#    fi
#    else
#    echo "ERROR($LINENO): Spaces not acceped in cell values"
#    fi
#else
#    echo "SKIP($LINENO): SQLITE3 not set"
#fi

if [[ -f "`which $PLINK`" ]]
then

    $PLINK \
        --make-bed \
        --bcf $BCF \
        --out $BCF.plink \
        --allow-extra-chr \
         2> /dev/null 1> /dev/null

    $PLINK \
         --bfile $BCF.plink \
         --freq \
         --allow-extra-chr \
         --out plink.out \
         2> /dev/null 1> /dev/null

    PLINK_COUNT=`cat plink.out.frq \
        | grep -v "CHR" \
        | awk '{if ($4 == "A") print $5*$6; else print (1-$5)*$6;}' \
        | awk '{s+=$1} END {print s}'`

    #GQT_COUNT=`$GQT query \
    GQT_COUNT=`$GQT query \
                    -i $BCF.gqt \
                    -d $DATA_PATH/more_fields.ped.db \
                    -p "" \
                    -g "maf()" \
                    | grep -v "^#" \
                    | awk '{
                        split($8,a,";");
                        for (i=1;i<=length(a);i++) {
                            split(a[i],b,"=");
                            if (b[1]=="GQT_0")
                                print b[2];
                        }
                    }' \
                | awk '{print $1*20}' \
                | awk '{s+=$1} END {print s}'`

    if [ $GQT_COUNT -eq $PLINK_COUNT ]
    then
        echo "SUCCESS($LINENO): GQT count matches PLINK count"
        rm  plink.out.frq plink.out.log plink.out.nosex
    else
        echo "ERROR($LINENO): GQT count does not match PLINK count. $GQT_COUNT vs $PLINK_COUNT"
        exit
    fi
else
    echo "SKIP($LINENO): PLINK not set"
fi

export BCFTOOLS_PLUGINS=":$HOME/src/bcftools/plugins/"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/src/htslib"
export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$HOME/src/htslib"

R=`bcftools plugin fill-AN-AC 2>&1`

if [[ $R == *"Could not load \"fill-AN-AC\""* ]]
then
    echo "BCFTOOLS plugin not found. Skipping remaining test. Please see README for help."
    exit
fi


BCFTOOLS_COUNT=`bcftools view $BCF\
    | bcftools plugin fill-AN-AC \
    | grep -v "#" \
    | awk '{
        split($8,a,";");
        for (i=1;i<=length(a);i++) {
            split(a[i],b,"=");
            if (b[1]=="AC")
                print b[2];
        }
    }' \
    | awk '{s+=$1} END {print s}'`

HET_COUNT=`$GQT query \
    -i $BCF.gqt \
    -d $DATA_PATH/more_fields.ped.db \
    -p "" \
    -g "count(HET)" \
    | grep -v "#" \
    | awk '{
        split($8,a,";");
        for (i=1;i<=length(a);i++) {
            split(a[i],b,"=");
            if (b[1]=="GQT_0")
                print b[2];
        }
    }' \
    | awk '{s+=$1} END {print s}'`

ALT_COUNT=`$GQT query \
    -i $BCF.gqt \
    -d $DATA_PATH/more_fields.ped.db \
    -p "" \
    -g "count(HOM_ALT)" \
    | grep -v "#" \
    | awk '{
        split($8,a,";");
        for (i=1;i<=length(a);i++) {
            split(a[i],b,"=");
            if (b[1]=="GQT_0")
                print b[2];
        }
    }' \
    | awk '{s+=$1*2} END {print s}'`

GQT_COUNT=`echo $HET_COUNT $ALT_COUNT | awk '{print $1+$2;}'`

if [ $GQT_COUNT -eq $BCFTOOLS_COUNT ]
then
    echo "SUCCESS($LINENO): GQT count matches BCFTOOLS count"
else
    echo "ERROR($LINENO): GQT count does not matche BCFTOOLS count. $GQT_COUNT vs $BCFTOOLS_COUNT"
    exit
fi

if [[ -f "`which $VCFTOOLS`" ]]
then
    if [[ -f "`which $BCFTOOLS`" ]]
    then
        BCFTOOLS view $BCF > $DATA_PATH/10.1e4.var.vcf
        echo -e "I0\nI1\nI2\nI3\nI4" > $DATA_PATH/A.txt
        echo -e "I5\nI6\nI7\nI8\nI9" > $DATA_PATH/B.txt
        $VCFTOOLS \
            --vcf $DATA_PATH/10.1e4.var.vcf \
            --weir-fst-pop $DATA_PATH/A.txt \
            --weir-fst-pop $DATA_PATH/B.txt \
            --out $DATA_PATH/A_vs_B \
        2> /dev/null > /dev/null
        tail -n+2 $DATA_PATH/A_vs_B.weir.fst | cut -f 3 > vcftools.fst.tmp
        $GQT fst \
            -i $BCF.gqt \
            -d $BCF.db \
            -p "BCF_Sample in ( 'I0', 'I1', 'I2', 'I3', 'I4')"\
            -p "BCF_Sample in ( 'I5', 'I6', 'I7', 'I8', 'I9')" \
        > $DATA_PATH/gqt.fst.vcf

        cat $DATA_PATH/gqt.fst.vcf \
        | grep -v "^#" \
        | awk '{
            split($8,a,";");
            for (i=1;i<=length(a);i++) {
                split(a[i],b,"=");
                if (b[1]=="GQT_fst")
                    print b[2];
            }
        }' \
        > gqt.fst.tmp

        MISSES=`paste vcftools.fst.tmp gqt.fst.tmp \
        | awk -F'\t' \
            'function abs(x){
                return ((x < 0.0) ? -x : x)
            } 
            {
                if (abs($1-$2) > 1e-05) print abs($0)
            }' \
        | wc -l`

        if [ $MISSES -eq 0 ]
        then
            echo "SUCCESS($LINENO): GQT Fst matches VCFTOOLS fst"
        else
            echo "ERROR($LINENO): GQT Fst does not match VCFTOOLS fst"
            cat $DATA_PATH/A_vs_B.weir.fst
            cat $DATA_PATH/gqt.fst.vcf
            exit
        fi
        rm $DATA_PATH/A_vs_B.weir.fst \
            $DATA_PATH/gqt.fst.vcf \
            gqt.fst.tmp \
            vcftools.fst.tmp
    fi
fi

clean_up
