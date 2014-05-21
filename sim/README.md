macs is a program that generates haplotypes.  To install:

    cd ~/src/
    wget https://macs.googlecode.com/files/macs-0.5d.tar.gz
    tar zxvf macs-0.5d.tar.gz
    cd macs
    make

macs creates haploids that need to be converted to diploids.  To do this we
generate 2N haploids then sum adjacent values to get N diploids:

    # simulate:
    # 100 individuals (200 haplotypes)
    # "genome" is 1Mb (1e6)
    # mutation and recombinaytion rate at 0.001
    ~/src/macs/macs 200 1e6 -T -t .001 -r .001 > 200.1e6.macs

    cat 200.1e6.macs | ~/src/genotq/sim/src/macs_hap_to_vari.py  > 100.1e6.vari

We can check the alt allele frequency of the variants by counting the number of
individuals that have that variant, then plot the histogram of those counts:

    ~/src/genotq/sim/src/get_alt_allele_freq.py \
        -d 100.1e6.vari \
        | ~/src/tools/plot/hist.py \
            -b 100 \
            -o 100_dip.png
