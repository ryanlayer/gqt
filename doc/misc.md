To convert all of the macs files to vari files
    ls *1e5*macs \
        | awk '
            {split($1,a,"."); 
            G=a[1]/2;
            print "cat "$1" | ~/src/genotq/sim/src/macs_hap_to_vari.py > "G"."a[2]".vari";}' \
        | sh &   

    ls *1e7*macs \
        | awk '
            {split($1,a,"."); 
            G=a[1]/2;
            print "cat "$1" | ~/src/genotq/sim/src/macs_hap_to_vari.py > "G"."a[2]".vari";}' \
        | sh &   

Convert vari files to indi files
    
    ls *vari \
        | sed -e "s/.vari//" \
        | awk '{ print "~/src/genotq/src/var_int_to_ind_int.py -v "$1".vari > "$1".indi &";}'

# 

Find the dimensions of the input file
    
    $ cd /Users/rl6sf/data/genotq/sim
    $ head -n1 10.1e4.vari | wc
        1      10      20
    $ wc -l 10.1e4.vari
        43 10.1e4.vari

Convert a varaint integer file to a variant uncompressed binary file

    $ ~/src/genotq/src/C/int_to_ubin \
        -i 10.1e4.vari \
        -f 10 \
        -r 43 \
        -o 10.1e4.varu

To convert an individual integer file just swap the ```-f``` and ```-r``` fileds

    $ head -n1 10.1e4.indi | wc
        1      43      86
    $ wc -l 10.1e4.indi
        10 10.1e4.indi
    $ ~/src/genotq/src/C/int_to_ubin \
        -i 10.1e4.indi \
        -f 43 \
        -r 10 \
        -o 10.1e4.indu

To get a genotype from an uncompressed binary files:

    $ ~/src/genotq/src/C/get_gt_from_ubin \
        -u 10.1e4.varu \
        -n 5 \
        -f 0,0,0,0,0 \
        -r 0,1,2,3,4
    0 0 2
    1 0 0
    2 0 1
    3 0 1
    4 0 0

Which matches the interger file:

    $ head -n 5 10.1e4.vari
    2 1 0 0 1 0 1 0 0 1
    0 0 0 0 0 0 0 0 1 0
    1 0 0 0 1 0 0 0 0 1
    1 0 0 0 1 0 0 0 0 0
    0 0 0 0 0 0 2 0 0 0

To OR a set of records 

    $ ~/src/genotq/src/C/or_gt_records_from_ubin \
        -u 10.1e4.varu \
        -n 3 \
        -r 0,1,2
    3 1 0 0 1 0 1 0 1 1

Which is correct

    $ head -n 3 10.1e4.vari
    2 1 0 0 1 0 1 0 0 1
    0 0 0 0 0 0 0 0 1 0
    1 0 0 0 1 0 0 0 0 1
 
    3 1 0 0 1 0 1 0 1 1


To OR a set of fields 

    $ ~/src/genotq/src/C/or_gt_fields_from_ubin \
        -u 10.1e4.indu \
        -n 3 \
        -f 0,1,2
    3 1 0 0 1 0 1 0 1 1

Which is correct

    $ cat 10.1e4.indi | cut -d" " -f1,2,3
    2 0 1
    1 0 0
    0 0 0
    0 0 0
    1 0 1
    0 0 0
    1 0 0
    0 0 0
    0 1 0
    1 0 1

