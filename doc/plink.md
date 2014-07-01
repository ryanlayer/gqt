# plink

##HapMap example

    cd ~/data/genotq/plink/hapmap_ex
    wget http://pngu.mgh.harvard.edu/~purcell/plink/hapmap1.zip
    unzip hapmap1.zip

    head hapmap1.map

    1 rs6681049 0 1
    1 rs4074137 0 2
    1 rs7540009 0 3
    1 rs1891905 0 4
    1 rs9729550 0 5
    1 rs3813196 0 6
    1 rs6704013 0 7
    1 rs307347 0 8
    1 rs9439440 0 9
    1 rs3128342 0 10 

The cols are: chrom, id, genetic distance, and bas-pair position

    head hapmap1.ped | cut -d " " -f1-20

    HCB181 1 0 0 1 1 2 2 2 2 2 2 1 2 2 2 2 2 2 2
    HCB182 1 0 0 1 1 2 2 1 2 2 2 1 2 1 2 2 2 2 2
    HCB183 1 0 0 1 2 2 2 1 2 2 2 1 2 1 1 2 2 2 2
    HCB184 1 0 0 1 1 2 2 1 2 2 2 1 1 2 2 2 2 2 2
    HCB185 1 0 0 1 1 2 2 1 2 2 2 2 2 2 2 2 2 2 2
    HCB186 1 0 0 1 1 2 2 2 2 2 2 1 1 2 2 2 2 2 2
    HCB187 1 0 0 1 1 2 2 2 2 2 2 1 2 1 2 2 2 2 2
    HCB188 1 0 0 1 1 2 2 1 2 2 2 1 1 2 2 2 2 2 2
    HCB189 1 0 0 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    HCB190 1 0 0 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2

The cols are: family id, individual id, paternal id, maternal id, sex,
phenotype

phenotype can be "missing" by setting it to -9

Genotypes (7+) are any character (e.g. 1,2,3,4 or A,C,G,T or anything else) except 0 which is, by default, the missing genotype character.

Making a binary PED file

    plink --file hapmap1 --make-bed --out hapmap1

### Transposed filesets:

transposed fileset, containing two text files: one (TPED)
containing SNP and genotype information where one row is a SNP; one (TFAM)
containing individual and family information, where one row is an individual.

You can generate transposed filesets with the --transpose option, described in
the data management section


### Extract a subset of individuals

    plink --file data --keep mylist.txt

where the file mylist.txt is just a list of Family ID / Individual ID pairs,
one set per line, i.e. one person per line. 

    head hapmap1.ped | cut -d " " -f1-2 > mylist.txt
    plink --bfile hapmap1 --freq --out freq_stat --keep mylist.txt
    head freq_stat.frq

    CHR       SNP   A1   A2          MAF  NCHROBS
    1   rs6681049    1    2            0       20
    1   rs4074137    1    2         0.25       20
    1   rs7540009    0    2            0       20
    1   rs1891905    1    2          0.5       20
    1   rs9729550    1    2          0.2       20
    1   rs3813196    1    2            0       20
    1   rs6704013    0    2            0       20
    1    rs307347    0    2            0       10
    1   rs9439440    0    2            0       20

    grep -n rs6681049 hapmap1.map

    1:1 rs6681049 0 1

    cat hapmap1.ped | head | cut -d " " -f7,8
    2 2
    2 2
    2 2
    2 2
    2 2
    2 2
    2 2
    2 2
    2 2
    2 2

    grep -n rs4074137 hapmap1.map
    2:1 rs4074137 0 2

    cat hapmap1.ped | head | cut -d " " -f9,10
    2 2
    1 2
    1 2
    1 2
    1 2
    2 2
    2 2
    1 2
    2 2
    2 2

## Convert plain text to plink-readable

We need both a map file and a ped file.  

### map file

The map file will contain 1 entry per variant.  The plain text file has the
number of records on the first line and the number of fields on the next, so
for a by-individual file (```*.ind.txt```), the number of variants will be on
the first line and by-variant (```*.ind.txt```) will be on the second.  

    cd ~/data/genotq/sim
    ~/src/genotq/src/py/plt_to_map.py \
        -p 100.1e8.ind.txt \
        -i \
        > 100.1e8.ind.map
    head 100.1e8.ind.map

    1 v0 0 1
    1 v1 0 2
    1 v2 0 3
    1 v3 0 4
    1 v4 0 5
    1 v5 0 6
    1 v6 0 7
    1 v7 0 8
    1 v8 0 9
    1 v9 0 10
    
### ped file

The ped file contains 1 line per individual.  The first 6 columns are family id,
individual id, paternal id, maternal id, sex, phenotype, and each individual
needs a unqiue family id + individual id.  The 7+ colums are pairs, where the
entries are the ids of the alleles.  The values can be 1,2,3,4 or A,C,T,G etc.
but we will just use 1 as ref and 2 as alt.  Here the mapping from plt to ped
is:
*   0 (homozygous ref)  -> 1 1
*   1 (heterozygous)  -> 1 2
*   2 (homozygous alt)  -> 2 2
*   3 (unknown )  -> 0 0
It only makes sense to use the by-individual plain text file here.

    cd ~/data/genotq/sim
    ~/src/genotq/src/py/ind_plt_to_ped.py \
        -p 100.1e8.ind.txt \
        > 100.1e8.ind.ped 
    head 100.1e8.ind.ped | cut -f1-30


    I1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1
    I2 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1
    I3 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1
    I4 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1
    I5 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 2 1 2 1 1 1 1 1 1
    I6 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1
    I7 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 2 1 2 1 1 1 1 1 1
    I8 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1
    I9 1 0 0 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1
    I10 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 2 1 2 1 2 1 1 1 1 1 1

    head -n 12 100.1e8.ind.txt | cut -d" " -f1-12

    588830
    100
    0 0 0 0 0 0 0 0 2 0 0 0
    0 0 0 0 0 0 0 0 1 0 0 0
    0 0 0 0 0 0 0 0 1 0 0 0
    0 0 0 0 0 0 0 0 1 0 0 0
    0 0 0 0 0 0 1 1 1 0 0 0
    0 0 0 0 0 0 0 0 2 0 0 0
    0 0 0 0 0 0 1 1 1 0 0 0
    0 0 0 0 0 0 0 0 1 0 0 0
    0 0 0 1 0 0 0 0 2 0 0 0
    0 0 0 0 0 1 1 1 1 0 0 0

    The partial first two lines are:
        
        0 0 0 0 0 0 0 0 2 0 0 0
        0 0 0 0 0 0 0 0 1 0 0 0

    Which should map to:
    
        1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1
        1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1

    Which matches the output above

### Test converted map and ped files

    plink --file 100.1e8.ind

    @----------------------------------------------------------@
    |        PLINK!       |     v1.07      |   10/Aug/2009     |
    |----------------------------------------------------------|
    |  (C) 2009 Shaun Purcell, GNU General Public License, v2  |
    |----------------------------------------------------------|
    |  For documentation, citation & bug-report instructions:  |
    |        http://pngu.mgh.harvard.edu/purcell/plink/        |
    @----------------------------------------------------------@

    Web-based version check ( --noweb to skip )
    Recent cached web-check found... OK, v1.07 is current

    +++ PLINK 1.9 is now available! See above website for details +++

    Writing this text to log file [ plink.log ]
    Analysis started: Tue Jun 17 12:30:55 2014

    Options in effect:
        --file 100.1e8.ind

    588830 (of 588830) markers to be included from [ 100.1e8.ind.map ]    
    100 individuals read from [ 100.1e8.ind.ped ]
    100 individuals with nonmissing phenotypes
    Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)
    Missing phenotype value is also -9
    0 cases, 100 controls and 0 missing
    100 males, 0 females, and 0 of unspecified sex
    Before frequency and genotyping pruning, there are 588830 SNPs
    100 founders and 0 non-founders found
    Total genotyping rate in remaining individuals is 1
    0 SNPs failed missingness test ( GENO > 1 )
    0 SNPs failed frequency test ( MAF < 0 )
    After frequency and genotyping pruning, there are 588830 SNPs
    After filtering, 0 cases, 100 controls and 0 missing
    After filtering, 100 males, 0 females, and 0 of unspecified sex

    Analysis finished: Tue Jun 17 12:31:51 2014

### Covert to binary format

    plink --file 100.1e8.ind --make-bed --out 100.1e8.ind

    @----------------------------------------------------------@
    |        PLINK!       |     v1.07      |   10/Aug/2009     |
    |----------------------------------------------------------|
    |  (C) 2009 Shaun Purcell, GNU General Public License, v2  |
    |----------------------------------------------------------|
    |  For documentation, citation & bug-report instructions:  |
    |        http://pngu.mgh.harvard.edu/purcell/plink/        |
    @----------------------------------------------------------@

    Web-based version check ( --noweb to skip )
    Recent cached web-check found... OK, v1.07 is current

    +++ PLINK 1.9 is now available! See above website for details +++

    Writing this text to log file [ 100.1e8.ind.log ]
    Analysis started: Tue Jun 17 12:33:09 2014

    Options in effect:
        --file 100.1e8.ind
        --make-bed
        --out 100.1e8.ind

    ** For gPLINK compatibility, do not use '.' in --out **
    588830 (of 588830) markers to be included from [ 100.1e8.ind.map ]    
    100 individuals read from [ 100.1e8.ind.ped ]
    100 individuals with nonmissing phenotypes
    Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)
    Missing phenotype value is also -9
    0 cases, 100 controls and 0 missing
    100 males, 0 females, and 0 of unspecified sex
    Before frequency and genotyping pruning, there are 588830 SNPs
    100 founders and 0 non-founders found
    Total genotyping rate in remaining individuals is 1
    0 SNPs failed missingness test ( GENO > 1 )
    0 SNPs failed frequency test ( MAF < 0 )
    After frequency and genotyping pruning, there are 588830 SNPs
    After filtering, 0 cases, 100 controls and 0 missing
    After filtering, 100 males, 0 females, and 0 of unspecified sex
    Writing pedigree information to [ 100.1e8.ind.fam ]
    Writing map (extended format) information to [ 100.1e8.ind.bim ]
    Writing genotype bitfile to [ 100.1e8.ind.bed ]
    Using (default) SNP-major mode

    Analysis finished: Tue Jun 17 12:34:07 2014

The file is for the compressed binary is almost identical to the genotq
compressed binary, both about 14M:

    -rw-r--r--+ 1 rl6sf  staff   14720808 Jun  3 14:31 100.1e8.ind.ubin
    -rw-r--r--+ 1 rl6sf  staff   14720753 Jun 17 12:34 100.1e8.ind.bed

At this size, the compressed binary is larger than the WAH (21M unsorted, 15M sorted):

    -rw-r--r--+ 1 rl6sf  staff   21762816 Jun  3 14:41 100.1e8.ind.wah
    -rw-r--r--+ 1 rl6sf  staff   15513828 Jun  3 16:37 100.1e8.ind.allele_sort.wah


### Runtimes for text/binary
    
    time plink --bfile 100.1e8.ind --freq --out freq_stat >/dev/null

    real    0m6.060s
    user    0m5.853s
    sys     0m0.198s

    time plink --file 100.1e8.ind --freq --out freq_stat >/dev/null

    real    0m58.628s
    user    0m58.229s
    sys     0m0.350s

### Try getting a subset

    head 100.1e8.ind.ped | cut -d " " -f1-2 > keep_list.txt
    time plink --bfile 100.1e8.ind --freq --out freq_stat --keep keep_list.txt >/dev/null 

    real    0m5.444s
    user    0m5.236s
    sys     0m0.203s
    

### Get all of the map files

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
