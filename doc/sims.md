
### 2014.04.23

    ssh radar
    cd /localtmp/rl6sf/genotq/sim

Create the vair files for a human-sized geno at at different population sizes
``SIZE=1e8``` and ```SIZE=3e9```

    POPS="100
    500
    1000
    5000
    10000"

    for POP in $POPS
    do
        HAP=`calc $POP*2`
        ~/src/macs/macs $HAP $SIZE -T -t .001 -r .001 2>/dev/null \
            | ~/src/genotq/sim/src/macs_hap_to_vari.py  \
            > $POP.$SIZE.vari &
    done

Convert vari to indi files:

    cd /localtmp/rl6sf/genotq/sim/small

    ~/src/genotq/src/py/var_int_to_ind_int.py -v 100.1e8.vari > 100.1e8.indi &
    ~/src/genotq/src/py/var_int_to_ind_int.py -v 500.1e8.vari > 500.1e8.indi &
    ~/src/genotq/src/py/var_int_to_ind_int.py -v 1000.1e8.vari > 1000.1e8.indi &
    ~/src/genotq/src/py/var_int_to_ind_int.py -v 5000.1e8.vari > 5000.1e8.indi &
    ~/src/genotq/src/py/var_int_to_ind_int.py -v 10000.1e8.vari > 10000.1e8.indi &
    
Covert integer files to uncompress binary files

    ~/src/genotq/scripts/vari2varu -i 100.1e8.vari -o 100.1e8.varu
    ~/src/genotq/scripts/vari2varu -i 500.1e8.vari -o 500.1e8.varu
    ~/src/genotq/scripts/vari2varu -i 1000.1e8.vari -o 1000.1e8.varu
    ~/src/genotq/scripts/vari2varu -i 5000.1e8.vari -o 5000.1e8.varu
    ~/src/genotq/scripts/vari2varu -i 10000.1e8.vari -o 10000.1e8.varu

    ~/src/genotq/scripts/vari2varu -i 100.1e8.indi -o 100.1e8.indu
    ~/src/genotq/scripts/vari2varu -i 500.1e8.indi -o 500.1e8.indu
    ~/src/genotq/scripts/vari2varu -i 1000.1e8.indi -o 1000.1e8.indu
    ~/src/genotq/scripts/vari2varu -i 5000.1e8.indi -o 5000.1e8.indu
    ~/src/genotq/scripts/vari2varu -i 10000.1e8.indi -o 10000.1e8.indu

Set some varibles:
    
    POPS="100
    500
    1000
    5000
    10000"

    GSIZE="1e8"

    TYPES="vari
    indi
    varu
    indu
    vari.gzip
    indi.gzip"


Gzip, bzip, and index (with grabix) all the text files:
    POPS="100
    500
    1000
    5000" 
    TXT_TYPES="vari
    indi"
    for POP in $POPS
    do
        for TXT_TYPE in $TXT_TYPES
        do
             gzip -c $POP.$GSIZE.$TXT_TYPE > $POP.$GSIZE.$TXT_TYPE.gz &
        done
    done

    for POP in $POPS
    do
        for TXT_TYPE in $TXT_TYPES
        do
             bgzip -c $POP.$GSIZE.$TXT_TYPE > $POP.$GSIZE.$TXT_TYPE.bgz &
        done
    done

    for POP in $POPS
    do
        for TXT_TYPE in $TXT_TYPES
        do
             grabix index $POP.$GSIZE.$TXT_TYPE.bgz &
        done
    done

Collection Metrics    

Size
    EXTS="vari
    indi
    vari.gz
    indi.gz
    vari.bgz
    indi.bgz
    varu
    indu"
 
    for EXT in $EXTS
    do
        for POP in $POPS
        do
            ls -l $POP.$GSIZE.$EXT | awk '{OFS="\t"; print $9,$5;}'
        done
    done

### 2014.05.09

    ssh radar
    cd /localtmp/rl6sf/genotq/sim/small

    ~/src/genotq/bin/or_gt_fields_from_ubin \
        -u 100.1e8.indu \
        -n 3 \
        -f 1,2,3
        
### 2014.05.12

    ssh radar
    cd /localtmp/rl6sf/genotq/sim/tiny
    RECORDS=`wc -l 10.1e4.indi`
    FIELDS=`head -n 1 10.1e4.indi  | awk '{print NF;}'`
    ~/src/genotq/bin/int_to_ubin \
        -i 10.1e4.indi \
        -r $RECORDS \
        -f $FIELDS \
        -o 10.1e4.indu

    ~/src/genotq/bin/or_gt_records_from_uint \
        -u 10.1e4.indi \
        -n 3 \
        -r 1,2,3 \
        -F 43 \
        -R 10
    ~/src/genotq/bin/or_gt_fields_from_uint \
        -u 10.1e4.vari \
        -n 3 \
        -f 1,2,3 \
        -F 43 \
        -R 10



    ~/src/genotq/bin/or_gt_records_from_ubin \
        -u 10.1e4.indu \
        -n 3 \
        -r 1,2,3 

    ~/src/genotq/bin/or_gt_fields_from_uint \
        -u 10.1e4.indi \
        -n 3 \
        -f 1,2,3 \
        -F 43 \
        -R 10

    ~/src/genotq/bin/or_gt_fields_from_ubin \
        -u 10.1e4.indu \
        -n 3 \
        -f 1,2,3 


    FILE=10000.1e8
    RECORDS=`cat $FILE.indi | wc -l`
    FIELDS=`head -n 1 $FILE.indi  | awk '{print NF;}'`

    ~/src/genotq/bin/int_to_ubin \
        -i $FILE.indi \
        -r $RECORDS \
        -f $FIELDS \
        -o $FILE.indu

    ~/src/genotq/bin/int_to_ubin \
        -i $FILE.vari \
        -r $RECORDS \
        -f $FIELDS \
        -o $FILE.varu



    bgzip -c $FILE.vari > $FILE.vari.bgz &
    bgzip -c $FILE.indi > $FILE.indi.bgz &
