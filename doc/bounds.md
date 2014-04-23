
Get simulated variants for 100 individuals with a 1e6 genome size

    $ cd ~/data/genotq/sim
    $ ~/src/macs/macs 200 1e6 -T -t .001 -r .001 > 200.1e6.macs
    $ cat 200.1e6.macs \
        | ~/src/genotq/sim/src/macs_hap_to_vari.py \
        > 100.1e6.vari

Verify that there are 100 genotypes

    $ head -n 1 100.1e6.vari | wc 
        1     100     200

Count the number of variants

    $ wc -l 100.1e6.vari
        5921 100.1e6.vari

And get the size

    $ ls -l 100.1e6.vari
    -rw-r--r--+ 1 rl6sf  staff  1184200 Apr  1 11:16 100.1e6.vari

Now convert the variable-based file (vari) to an individual-based file (indi)

    $ ~/src/genotq/src/var_int_to_ind_int.py -v 100.1e6.vari > 100.1e6.indi
    $ wc -l 100.1e6.indi
        100 100.1e6.indi
    $ head -n 1 100.1e6.indi | wc 
        1    5921   11842

Compress and check the sizes

    $ gzip -9 -c 100.1e6.indi > 100.1e6.indi.gz
    $ gzip -9 -c 100.1e6.vari > 100.1e6.vari.gz
    $ ls -l 100.1e6.vari 100.1e6.vari.gz 100.1e6.indi 100.1e6.indi.gz
        -rw-r--r--+ 1 rl6sf  staff  1184200 Apr 15 12:52 100.1e6.indi
        -rw-r--r--+ 1 rl6sf  staff   100732 Apr 15 12:53 100.1e6.indi.gz
        -rw-r--r--+ 1 rl6sf  staff  1184200 Apr  1 11:16 100.1e6.vari
        -rw-r--r--+ 1 rl6sf  staff    59990 Apr 15 12:55 100.1e6.vari.gz

Compression ratio of indi/indi.gz

    $ calc 100732/1184200
    0.0850633338963013

Compression ratio of vari/vari.gz

    $ calc 59990/1184200
    0.0506586725215335

Ratio of indi.gz/vari.gz
    
    $ calc 100732/59990
    1.67914652442074    

vari compresses better than indi by a factor of 1.67

Now we need to test the runtime of queiries against the two types.

Query: select variants that are non-reference among a given set of individuals.

Test to see if the program works:

    $ cat 100.1e6.indi | \
        ~/src/genotq/src/get_common_genot.py \
            -i -I 1,2,3,4 \
        | summation
    150885

    $ cat 100.1e6.vari | \
        ~/src/genotq/src/get_common_genot.py \
            -v -I 1,2,3,4 \
        | summation
    150885


Now we need to time the two:

    $ time cat 100.1e6.vari \
        | ~/src/genotq/src/get_common_genot.py \
            -v -I 1,2,3,4,5,6,7,8,9 >/dev/null
    real    0m0.181s
    user    0m0.136s
    sys     0m0.046s

    $ time cat 100.1e6.indi \
        | ~/src/genotq/src/get_common_genot.py \
            -i -I 1,2,3,4,5,6,7,8,9 >/dev/null
    real    0m0.157s
    user    0m0.111s
    sys     0m0.046s

To get a better estimate of time, a larger test is needed.  For a given 
individual/variant integer pair of genotype files, we find how long it takes to
get the common alternate alleles among a randomly selected number of people.
In this case we get 10 people 10 times, and return the mean.

    ~/src/genotq/test/speed_get_common_genot_py.sh \
        -i 100.1e6.indi \
        -v 100.1e6.vari \
        -N 10 \
        -n 10 \
        -t 100
    100.1e6.vari 0.1805,0.00335410196624969
    100.1e6.indi 0.1606,0.00300665927567458

Try a larger test

    ~/src/genotq/test/speed_get_common_genot_py.sh \
        -i 1000.1e6.indi \
        -v 1000.1e6.vari \
        -N 10 \
        -n 10 \
        -t 100
    1000.1e6.vari 0.357,0.00538516480713451
    1000.1e6.indi 0.1799,0.00284429253066558

Scan the sizes
    
     ~/src/genotq/test/scan_speed_get_common_genot_py.sh \
        -i 1000.1e6.indi \
        -v 1000.1e6.vari \
        -N 3 \
        -S 10 \
        -E 100 \
        -X 10 \
        -t 100
    
    10 1000.1e6.vari 0.359333333333333,0.0060184900284226 1000.1e6.indi 0.176666666666667,0.000471404520791032
    20 1000.1e6.vari 0.412,0.00355902608401044 1000.1e6.indi 0.249333333333333,0.00478423336480245
    30 1000.1e6.vari 0.481,0.00355902608401044 1000.1e6.indi 0.314,0.00432049379893858
    40 1000.1e6.vari 0.541333333333333,0.0103387082795139 1000.1e6.indi 0.171333333333333,0.00188561808316412
    50 1000.1e6.vari 0.594666666666667,0.00679869268479039 1000.1e6.indi 0.164333333333333,0.0016996731711976
    60 1000.1e6.vari 0.657666666666667,0.00590668171555646 1000.1e6.indi 0.157333333333333,0.00659966329107445
    70 1000.1e6.vari 0.726333333333333,0.00329983164553722 1000.1e6.indi 0.123666666666667,0.000471404520791032
    80 1000.1e6.vari 0.783,0.00509901951359279 1000.1e6.indi 0.131,0
    90 1000.1e6.vari 0.849,0.0155563491861041 1000.1e6.indi 0.130666666666667,0.000942809041582064
    100 1000.1e6.vari 0.907666666666667,0.0108423039781937 1000.1e6.indi 0.133,0.00282842712474619

     ~/src/genotq/test/scan_speed_get_common_genot_py.sh \
        -i 1000.1e6.indi \
        -v 1000.1e6.vari \
        -N 3 \
        -S 10 \
        -E 1000 \
        -X 10 \
        -t 100 \
        > get_common_genot_py_S10_E1000_X10_1000.1e6.indi_1000.1e6.vari.times
 
     ~/src/genotq/test/scan_speed_get_common_genot_c.sh \
        -i 1000.1e6.indi \
        -v 1000.1e6.vari \
        -N 3 \
        -S 10 \
        -E 1000 \
        -X 10 \
        -t 100 \
        > get_common_genot_c_S10_E1000_X10_1000.1e6.indi_1000.1e6.vari.times
 
     ~/src/genotq/test/scan_speed_get_common_genot_c.sh \
        -i 1000.1e5.indi \
        -v 1000.1e5.vari \
        -N 3 \
        -S 10 \
        -E 1000 \
        -X 10 \
        -t 100 \
        > get_common_genot_c_S10_E1000_X10_1000.1e5.indi_1000.1e6.vari.times
 
     ~/src/genotq/test/scan_speed_get_common_genot_c.sh \
        -i 1000.1e7.indi \
        -v 1000.1e7.vari \
        -N 3 \
        -S 10 \
        -E 1000 \
        -X 10 \
        -t 100 \
        > get_common_genot_c_S10_E1000_X10_1000.1e7.indi_1000.1e6.vari.times
 

Graph the results

    ~/src/genotq/graph/scan_speed_get_common_genot.py \
        -d ~/data/genotq/sim/get_common_genot_S10_E1000_X10_1000.1e6.indi_1000.1e6.vari.times \
        -o get_common_genot_S10_E1000_X10_1000.1e6.indi_1000.1e6.vari.times.png \
        --y_max 10 \
        --y_min 1 \
        --x_max 1000 \
        --x_min 1 \
        --minor 100:1600:200
   
