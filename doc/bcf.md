# VCF

### To convert from by-variant plain text to VCF:

    cd ~/data/genotq/sim
    ~/src/genotq/src/py/var_plt_to_vcf.py \
        -v 100.1e8.var.txt \
        > 100.1e8.var.vcf


### Then to convert to BCF

    bcftools view -O b 100.1e8.var.vcf > 100.1e8.var.bcf

### Compare file sizes

    -rw-r--r--+ 1 rl6sf  staff   117766011 Jun 17 13:37 100.1e8.var.txt
    -rw-r--r--+ 1 rl6sf  staff   257096641 Jun 17 13:43 100.1e8.var.vcf
    -rw-r--r--+ 1 rl6sf  staff    11550866 Jun 17 13:44 100.1e8.var.bcf

### Get all of the bcf files

    for f in `ls *ind.txt`
    do
        p=`echo $f|cut -d"." -f1-3`
        echo $p
        ~/src/genotq/src/py/plt_to_map.py \
            -p $p.txt \
            -i \
        > $p.map
    done

### Get all of the ped files

    for f in `ls *ind.txt`
    do
        p=`echo $f|cut -d"." -f1-3`
        echo $p
        ~/src/genotq/src/py/ind_plt_to_ped.py \
            -p $p.txt \
        > $p.ped 
    done

### Get all of the bed files

    for f in `ls *ind.ped`
    do
        p=`echo $f|cut -d"." -f1-3`
        echo $p
        plink --file $p --make-bed --out $p
    done
