To convert from VCF to a plain text files 

    zcat ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz \
        | ~/src/genotq/bin/vcf_to_plt \
            -f 1092 \
            -r 494328 \
        > chr22.var.txt


