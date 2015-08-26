GQT
========
Genotype Query Tools (GQT) is command line software and a C API for indexing 
and querying large-scale genotype data sets like those produced by 1000 Genomes,
the UK100K, and forthcoming datasets involving millions of genomes. GQT 
represents genotypes as compressed bitmap indices, which reduce 
computational burden of variant queries based on sample genotypes, 
phenotypes, and relationships by orders of magnitude over standard "variant-centric"
indexing strategies. This index can significantly expand the capabilities of 
population-scale analyses by providing interactive-speed queries to data sets 
with millions of individuals. 

[![Install and Demo Video](http://layerlab.org/gqt/GQT_SS.png)](http://www.youtube.com/watch?v=floxSA2OoM8)


## Table of Contents
1. [Installation] (#installation)
2. [Example work flow] (#example-workflow)
2. [Extra detail] (#extra-detail)

## Installation
GQT depends on htslib, sqlite3, and lex (flex).

*Step 1*. Install htslib.
```
git clone https://github.com/samtools/htslib.git
cd htslib
make
cd ..
```

*Step 2*. Download sqlite amalgamation source.
```
wget http://www.sqlite.org/2014/sqlite-amalgamation-3080701.zip
unzip sqlite-amalgamation-3080701.zip
```

*Step 3*. Check to see if your system has lex/flex installed.  If not, install.
```
lex -V
flex -V
```

If both fail then install from source, otherwise skip to step 4 
```
wget http://downloads.sourceforge.net/project/flex/flex-2.5.39.tar.bz2
bunzip2 flex-2.5.39.tar.bz2
tar xvf flex-2.5.39.tar
cd flex-2.5.39
./configure
make
make install
cd ..
```

*Step 4*. Get GQT source then modify the GQT Makefile by setting the
`HTS_ROOT` and `SQLITE_ROOT` variable in `src/Makfile` to reflect their
locations based on the directories in which they were placed during steps 1 and 2.  Compile.

```
git clone https://github.com/ryanlayer/gqt.git
cd gqt/
make
```

At this point, it is recommended that you copy the gqt binary to a directory that is on your PATH.



*Step 5 (Optional)* Run the GQT functional tests.

5a. Install bcftools (not necessary for GQT to function, but useful for the functional tests below).
```
cd ..
git clone https://github.com/samtools/bcftools
cd bcftools
make
```

5b. Install plink (v1.9) (not necessary for GQT to function, but useful for the functional tests below)..
```
# Download the appropriate binary from:
https://www.cog-genomics.org/plink2
# Now copy the plink binary to a directory on your PATH
```

5c. In addition, after you install [bcftools](https://github.com/samtools/bcftools),
you need to also update the directory assigned to [BCFTOOLS_PLUGIN](https://github.com/ryanlayer/gqt/blob/master/test/func/functional_tests.sh#L510) in your `gqt/test/func/functional_tests.sh`
file to be the plugins directory in the bcftools source tree. For example, if you compiled bcftool
in `~/src`, the correct path for BCFTOOLS_PLUGIN would be `~/src/bcftools/plugins`. That is:

    export BCFTOOLS_PLUGINS="$HOME/src/bcftools/plugins"

5d. Lastly, you may need to also upate the directories assigned to either [LD_LIBRARY_PATH](https://github.com/ryanlayer/gqt/blob/master/test/func/functional_tests.sh#L511) (linux)
or [DYLD_LIBRARY_PATH](https://github.com/ryanlayer/gqt/blob/master/test/func/functional_tests.sh#L512) (Mac OS X) in your `gqt/test/func/functional_tests.sh` file to include the directory where `libhts.so.1` is located. For example, if you compiled htslib in `~/src`, the correct settings would be:
    
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/src/htslib"
    export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$HOME/src/htslib"

Okay, now you are ready to run the tests. Make sure you are back in the make gqt source directory.

```
cd ../gqt/test/func
./functional_tests.sh
cd ../../..
```

## Example workflow

#### Create a GQT index from a BCF
GQT starts with a target VCF, VCF.GZ, or BCF file.  Download an example BCF
file.
```
wget --trust-server-names http://bit.ly/gqt_bcf
```

To index the BCF, GQT needs to know the number of columns and rows in the BCF.
This information can be either supplied with command line parameters or GQT can
extract this info from the CSI index produced by bcftools.
```
bcftools index chr11.11q14.3.bcf
```

With the BCF index created, now create the GQT index.
```
gqt/bin/gqt convert bcf -i chr11.11q14.3.bcf
```

This will create three files:
- `chr11.11q14.3.bcf.gqt` is the GQT index
- `chr11.11q14.3.bcf.bim` is the compressed variant metadata
- `chr11.11q14.3.bcf.vid` is the variant order information

#### Create a GQT sample database
GQT requires a database that describes the samples in the target BCF file.
This database is based on the samples in the BCF, and can be optionally
augmented by file (typically called PED file) that contains tab separated
fields, one line that gives the field name, and the subsequent lines give the
field values for each sample.  One column in this file must correspond to a
sample in the BCF file.  By default column 2 is used, but this can be modified
by using the `-c` option.  A warning will be generated for each record in the
PED file that does not have a corresponding sample in the VCF/BCF.  NOTE: Each
field in the first line will become a database column.  Given the naming
restriction, many special characters are not allowed (e.g., %, $, :), and all
spaces will be converted to underscores ("_").

For example, the following file includes Family ID, Population and Gender
information.  Using the default sample column (2), the each Individual ID will
be matched with the Sample ID in the VCF/BCF.

```
Family ID   Individual ID   Population  Gender
HG00096     HG00096         GBR         1
HG00171     HG00171         FIN         2
```

Without a PED file the database can be created from the BCF file:
```
gqt/bin/gqt convert ped -i chr11.11q14.3.bcf
```
This creates a database named `chr11.11q14.3.bcf.db`

When a PED file is available it can be included, thus allowing more 
sophisticated queries based on the phenotypes, ancestries, and 
relationships of the samples. More on this below.
```
wget --trust-server-names http://bit.ly/gqt_ped
gqt/bin/gqt convert ped -i chr11.11q14.3.bcf -p 1kg.phase3.ped
```
This creates a database named `1kg.phase3.ped.db`

#### Query the GQT index
When the database is based only on the BCF file (`chr11.11q14.3.bcf.db`),
queries can be based on the sample names. For example, you can find the
variants where at least two of the three listed individuals is heterozygous.
```
gqt/bin/gqt query \
    -i chr11.11q14.3.bcf.gqt \
    -p "BCF_Sample in ('NA21126','NA21127','NA21128')" \
    -g "count(HET)>1" \
    > NA21126-8.HET.vcf
```

When the database also includes the PED file (`1kg.phase3.ped.db`), it is
possible to query based on fields like gender and population.  For example, the
following query will find all variants that has a frequency of more than 10%
among the GBR population:
```
gqt/bin/gqt query \
    -i chr11.11q14.3.bcf.gqt \
    -d 1kg.phase3.ped.db \
    -p "Gender = 1 and Population ='GBR'" \
    -g "maf()>0.1" \
    > GBR.vcf
```

GQT can support any number of queries, where each further filters the variants
that are returned.  This is particular helpful for case control studies,
pedigree analysis, and for comparing populations.
```
gqt/bin/gqt query \
    -i chr11.11q14.3.bcf.gqt \
    -d 1kg.phase3.ped.db \
    -p "Gender = 1 and Population ='GBR'" \
    -g "maf()>0.1" \
    -p "Population = 'YRI'" \
    -g "maf()<0.1"
```

If you are only interested in how many variants pass a give set of filters,
then  the `-c` option can be used.  This is helpful for exploration because is
does not create a VCF file and therefore tends to be orders of magnitude
faster.
```
gqt/bin/gqt query \
    -i chr11.11q14.3.bcf.gqt \
    -d 1kg.phase3.ped.db \
    -p "Gender = 1 and Population ='GBR'" \
    -g "maf()>0.1" \
    -c
```

#### Integration with other variant tools
Since the output is VCF we can further filter variants by total depth with
bcftools:
```
gqt/bin/gqt query \
    -i chr11.11q14.3.bcf.gqt \
    -d 1kg.phase3.ped.db \
    -p "Gender = 1 and Population ='GBR'" \
    -g "maf()>0.1" \
    | bcftools view - -i 'DP>10000'
```

Or use bedtools to intersect the output with other genome annotation:
```
wget --trust-server-names http://bit.ly/gqt_genes
gqt/bin/gqt query \
    -i chr11.11q14.3.bcf.gqt \
    -d 1kg.phase3.ped.db \
    -p "Gender = 1 and Population ='GBR'" \
    -g "maf()>0.1" \
    | bedtools intersect -a genes.bed -b stdin
```

## Extra detail
GQT takes a BCF as input and produces three files, a (very small) compressed
index (.gqt), a compressed summary of the variant metadata (.bim), and a
variant ID file (.vid). The `convert` process rotates the data, sorts it by
alternate allele frequency, converts it to a bitmap index, then finally
compresses the data using the Word Aligned Hybrid (WAH) encoding
(http://dl.acm.org/citation.cfm?id=1132864).  Rotating the typical
variant-row/individual-column form of a BCF such that individuals become rows
allows for more efficient memory access patterns for individual based queries
(e.g., find the alternate allele frequency count for some set of individuals
within the cohort).  Sorting columns (variants) improves compression.
Converting to a bitmap index allows queries to be resolved using bit-wise
logical operations, which can compare 32 genotypes in a single fast operation.
WAH encoding gives near-optimal compression while allowing bit-wise logical
operations without inflation.  The combined result of these steps creates an
index that is only a fraction of the source BCF size, and queries against that
index complete in seconds.

The following [slides](http://quinlanlab.org/pdf/presentations/gtqGI2014v6.pdf)
provides a conceptual overview, as well as more details about the choice of
bitmaps and the word-aligned hybrid strategy for this problem.


### GQT Query Interface
#### Phenotype queries
Phenotype queries (specified with a `-p` option) are be resolved using data
from the sample database, and must contain expressions that are based on the
fields in that database.  When the database was built with `1kg.phase3.ped`
which includes files like Gender and Population it is possible to query for:

* Women:
```
-p "Gender==2"
```
* African Caribbean in Barbados:
```
-p "Population=='ACB'"
```
* Men of African decent: 
```
-p "Gender==1 and Population in ('YRI','LWK','GWD','MSL','ESN','ASW','ACB')"
```
* Men NOT of European decent: 
```
-p "Gender==1 and Population not in ('CEU', 'EUR', 'TSI', 'FIN', 'GBR', 'IBS')
```

#### Genotype queries
Genotype queries (specified with a `-g` option) are resolved using the GQT
index, and are composed of a set of predefined values, functions, and binary
relations:

* values
  * `HOM_REF`
  * `HET`
  * `HOM_ALT`
* functions
  * `count(HET)` The number of samples that are heterozygous
  * `pct(HET HOM_ALT)` The percent of samples that are either heterozygous or
    homozygous alternate 
  * `maf()` The minor allele frequency
* and binary relations 
  * `count(HET) => 1` at least one sample is heterozygous
  * `pct(HET HOM_ALT) < 0.1` less that 10% of individuals are either
    heterozygous or homozygous alternate

#### Output
The default output of GQT queries is VCF.  An extra line is added to the header
and to each line that gives the results of the query.  For example, a query
that include two queries
```
-p "Population=='GBR'" -g "count(HET HOM_ALT)" 
-p "Population=='FIN'" -g "pct(HET HOM_ALT)"
```

Then the header will include two results, one for the count and one for the percent:

```
##INFO=<ID=GTQ_0,Number=1,Type=Integer,Description="GQT count result from query 0">
##INFO=<ID=GTQ_1,Number=1,Type=Float,Description="GQT percent result from query 1">
```

And each line will have two corresponding values in the INFO field.  For example
`GTQ_0=2;GTQ_1=0.010753` indication that 2 individuals from the GBR
population and 0.01075 percent of the FIN population had a non-ref allele for
that variant.
