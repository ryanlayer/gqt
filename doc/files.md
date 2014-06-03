# Plain text file (plt)

# Uncompressed binary file

# Bitmap index

# WAH compressed bitmap index

# VCF

    

# BCF

# Grabix

# plink

HapMap example

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
