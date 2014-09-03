# 1000 Genomes Test

Get the data

    ssh radar
    cd /localtmp/rl6sf/genotq/1kg
    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
    gunzip ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz

Find the number of fields

    grep -m 1 "^#CHROM" ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf | cut -f 10- | awk '{print NF;}'  

    1092

Find the number of records

    grep -v "^#" ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf | wc -l

    855166

Convert to plt

    gtq convert vcf-plt \
        -i ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf \
        -o ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.var.txt \
        -r 855166 \
        -f 1092

Invert to ubin

    gtq convert plt-invert-ubin \
        -i ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.var.txt \
        -o ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.ind.ubin

Back to plt then sort

    gtq convert ubin-plt \
        -i ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.ind.ubin \
        -o ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.ind.txt
    gtq sort plt-field-freq \
        -i ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.ind.txt \
        -o ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.allele_sort.ind.txt

Back to ubin then wahbm

    gtq convert plt-ubin \
        -i ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.allele_sort.ind.txt \
        -o ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.allele_sort.ind.ubin

    gtq convert ubin-wahbm \
        -i ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.allele_sort.ind.ubin \
        -o ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.allele_sort.ind.wahbm

BGZip it

    bgzip \
        -c ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.ind.txt \
        > ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.ind.txt.bgz


Query it:


    N=1000
    R=`seq 1 999  | tr '\n' ','`
    R=`echo $R\1000`
    PRE="ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.allele_sort.ind"

Test for correctness

    gtq count plt -o gt -q 0 -i $PRE.txt -n $N -r $R | tr ' ' '\n' | summation
    7893594

    gtq count ubin -o gt -q 0 -i $PRE.ubin -n $N -r $R | tr ' ' '\n' | summation
    7893594

    gtq count wahbm -o gt -q 0 -i $PRE.wahbm -n $N -r $R | tr ' ' '\n' | summation
    7893594
    
Time it:


    gtq count plt -o gt -q 0 -i $PRE.txt -n $N -r $R -Q -t
    3850470

    gtq count ubin -o gt -q 0 -i $PRE.ubin -n $N -r $R -Q -t
    4160921

    gtq count wahbm -o gt -q 0 -i $PRE.wahbm -n $N -r $R -Q -t
    3191647


Do the same process for chr15, which is pretty close to the size of our simulation

    grep -m 1 "^#CHROM" ALL.chr15.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf | cut -f 10- | awk '{print NF;}'
    1092

    grep -v "^#" ALL.chr15.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf | wc -l
    1130554

    gtq convert vcf-plt -i ALL.chr15.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf -o 1kg.chr15.var.txt -f 1092 -r 1130554
    gtq convert plt-vcf -i 1kg.chr15.var.txt -o 1kg.chr15.var.vcf
    gtq convert plt-invert-ubin -i  1kg.chr15.var.txt -o 1kg.chr15.ind.ubinI
    head -n 2 1kg.chr15.var.txt
      1092
      1130554
    echo -e "1130554\n1092" > 1kg.chr15.ind.txt
    gtq view ubin -i 1kg.chr15.ind.ubin >> 1kg.chr15.ind.txt
    gtq sort plt-field-freq -i 1kg.chr15.ind.txt -o 1kg.chr15.ind.allele_sort.txt
    gtq convert plt-ubin -i 1kg.chr15.ind.allele_sort.txt -o 1kg.chr15.ind.allele_sort.ubin
    gtq convert ubin-wahbm -i 1kg.chr15.ind.allele_sort.ubin -o 1kg.chr15.ind.allele_sort.wahbm

    bcftools view -O b 1kg.chr15.var.vcf -o 1kg.chr15.var.bcf

    ~/src/genotq/src/py/ind_plt_to_ped.py -p 1kg.chr15.ind.allele_sort.txt > 1kg.chr15.ind.allele_sort.ped
    ~/src/genotq/src/py/plt_to_map.py -i -p 1kg.chr15.ind.allele_sort.txt > 1kg.chr15.ind.allele_sort.map
    ~/src/plink2/plink --file 1kg.chr15.ind.allele_sort --make-bed --out 1kg.chr15.ind.allele_sort
     
