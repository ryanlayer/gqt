[toc]

Overview
========

NOTE: GQT is in an alpha state. 

Genome Query Tools (GQT) is command line software and a C API for storing 
and querying large-scale genotype data sets like those produced by 1000 Genomes,
the Uk100K, and forthcoming datasets involving millions of genomes. GQT 
represents genotypes as compressed bitmap indices, which reduce the storage and
compututational burden by orders of magnitude. This index can significantly expand 
the capabilities of population-scale analyses by providing interactive-speed
queries to data sets with millions of individuals. An example work flow:

Use GQT to create a sample-based index of a BCF file containing 1092
individuals and 494328 variants.

    gqt convert bcf \
        -r 494328 \
        -f 1092 \
        -i 1kg.chr22.bcf 

BCF is preferred for speed, but the same command works for VCF & compressed VCF:

    gqt convert bcf \
        -r 494328 \
        -f 1092 \
        -i 1kg.chr22.vcf 

    gqt convert bcf \
        -r 494328 \
        -f 1092 \
        -i 1kg.chr22.vcf.gz

Now, create a database of the PED file describing the phenotypes and relationships
of the 1092 samples. This enables powerful queries based on sample phenotype status,
ancestry, relationships, or any attribute you want, so long as it is the PED file.
NOTES: 1) The file doesn't strictly need to be PED, but that is a rational format 
to use. 2) The number and order of the individuals in the "PED" must exactly match
the number and order of the individuals in the VCF/VFC.GZ/BCF file.

    gqt convert ped \
        -i 1kg.ped 

Now that we have a GQT index and a sample database, we can query the index
for variants that meet query criteria based on both sample genotype and 
phenotype. For example, the following command will find all variants where 
at least 10 GBR individuals are heterozygous:

    gqt query \
        -i 1kg.chr22.bcf.gqt \
        -d 1kg.ped.db \
        -p "Population = 'GBR'" \
        -g "count(HET) >= 10"

Or, find all variants where the alternate allele frequency >= 10% affected 
in GBR individuals.

    gqt query \
        -i 1kg.chr22.bcf.gqt \
        -d 1kg.ped.db \
        -p "Population = 'GBR'" \
        -g "freq(HET) >= 0.10"

Or, find all variants where _all_ affected individuals are homozygous.

    gqt query \
        -i 1kg.chr22.bcf.gqt \
        -d 1kg.ped.db \
        -p "Phenotype = 2" \
        -g "HOMO_REF”

The default output of GQT queries is VCF with genotypes. As such, we can 
further filter variants by total depth with bcftools:

    gqt query \
        -i 1kg.chr22.bcf.gqt \
        -d 1kg.ped.db \
        -p "Population = 'GBR'" \
        -g "count(HET) >= 10" \
    | bcftools view - -i 'DP>50000'

A bit more detail.
=====================
GQT takes a BCF as input and produces three files, a (very small) compressed
index (.gqt), a compressed summary of the variant metadata (.bim), and a
variant ID file (.vid). The `convert` process rotates the data, sorts it by alternate
allele frequency, converts it to a bitmap index, then finally compresses the
data using the Word Aligned Hybrid (WAH) encoding
(http://dl.acm.org/citation.cfm?id=1132864).  Rotating the typical
variant-row/individual-column form of a BCF such that individuals become rows
allows for more efficient memory access patterns for individual based queries
(e.g., find the alternate allele frequency count for some set of individuals
within the cohort).  Sorting columns (variants) improves compression.
Converting to a bitmap index allows queries to be resolved using bit-wise
logical operations, which can compare 32 genotypes in a single fast operation.
WAH encoding gives near-optimal compression while allowing bit-wise logical
operations without inflation.

The combined result of these steps creates an index that is only a fraction of
the source BCF size, and queries against that index complete in seconds.  In
the case of chr22 from 1000 Genomes, the BCF file is 11G and the WAH index is
42M, and the alternate allele frequency count for 100 of those individuals can
be found in 0.188 seconds.

The following [slides](http://quinlanlab.org/pdf/presentations/gtqGI2014v6.pdf)
provides a conceptual overview, as well as more details about the choice of
bitmaps and the word-aligned hybrid strategy for this problem.

GQT queries are individual-focused.  That is, results such as alternate allele
frequency count or the set of non-reference alleles are based on a subset of
individuals from the population.  Identifying the target individuals can either
be done manually by selecting the individual’s zero-based ID through the
command line argument `-r`, or by querying an associated PED file
(http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped).  To query a PED
file, you must first convert the PED file to an sqlite3
(http://www.sqlite.org/) database by running the `gqt convert ped` commad
(see below).  Once this database has been created, simply pass the conditions
that define your target population population (e.g., `Population = 'GBR' AND
Paternal ID ='NA19679'`).  NOTE:  all white space characters in the PED file
will be converted to an underscore in the database (e.g., `Paternal ID` becomes
`Paternal_ID`), and all queries must follow this syntax.

Installation
============
GQT depends on htslib, sqlit3, and lex (flex).

*Step 1*. Install htslib.

    $ git clone https://github.com/samtools/htslib.git
    $ cd htslib
    $ make

*Step 2*. Download sqlite amalgamation source.

    $ wget http://www.sqlite.org/2014/sqlite-amalgamation-3080701.zip
    $ unzip sqlite-amalgamation-3080701.zip

*Step 3*. Check to see if your system has lex/flex installed.  If not, install.

    $ lex -V
    $ flex -V
    # if both fail then install, otherwise skip to step 4 
    $ wget http://downloads.sourceforge.net/project/flex/flex-2.5.39.tar.bz2
    $ bunzip2 flex-2.5.39.tar.bz2
    $ tar xvf flex-2.5.39.tar
    $ cd flex-2.5.39
    $ ./configure
    $ make
    $ make install

*Step 4*. Modify the GQT Makefile by setting the `HTS_ROOT` and `SQLITE_ROOT`
variable in `src/c/Makfile` to reflect their locations.

*Step 5*. Compile GQT

    $ cd gqt/
    $ make

*Step 6*. Test GQT.  First, make the same changes to `HTS_ROOT` and
`SQLITE_ROOT` in `test/unit/Makefile`. To test AVX2 instruction, uncomment the
`AVX` variable definition.

    $ cd test/unit
    $ make
    $ cd ../func
    $ bash functional_tests.sh


Example usage with 1000 Genomes data.
======================================

In the example below, we demonstrate how to create a GQT index of a compressed
VCF file from the 1000 Genomes Project. We then demonstrate how to use this
index to quickly query variants based on sample genotypes. Please note that the
conversion of compressed VCF to a GQT index will take only 4 or 5 minutes.
However, if you are impatient, we have posted all of the converted files
described above on our lab
[website](http://quinlanlab.cs.virginia.edu/gqt-example/).

*Step 1*. Download the Phase 1 1000 genomes chr22 file, and move the giant file
name to something more manageable.

	$ wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
        $ mv ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz 1kg.chr22.vcf.gz

*Step 2*. Use GQT to make an index of the genotypes in the compressed VCF file
using the `gqt covert bcf` command.  While this command is labeled `bcf`, it can
also accept `vcf` and `vcf.gz`.  Currently, one must tell GQT how many variants
(`-r`) and how many samples (`-f`) are in the file. Future versions will detect
this automatically.  To extract this info, one can use the
[bcftools](http://samtools.github.io/bcftools/) `stats` command. 

	# use bcftools stats to get the number of variants and individuals in the file.
	$ bcftools stats 1kg.chr22.vcf.gz | \
	    grep SN | \
	    head -4
	# SN, Summary numbers:
	# SN	[2]id	[3]key	[4]value
	SN	0	number of samples:	1092
	SN	0	number of records:	494328

As you can see, there are 494328 variants (records) and 1092 samples. We
provide this information to the GQT `convert` tool and ask it to make an index
of the compressed VCF file. In addition, the tool will create a `.bim` file
that stores information about the variants in the original file, and a `.vid`
file that stores the record number of the variants. GQT can very quickly report
results in valid VCF information from the `.bim`.  These results will not
included the sample columns.  GQT can retrieved that extra data be reaching
back into the source BCF (or VCF) file using the `.vid` file.  We will show
examples of both cases below.

        $ gqt convert bcf        \
	    -r 494328             \
	    -f 1092               \
	    -i 1kg.chr22.vcf.gz

*Step 3*.  Query the GQT index.  In this case we are going to simply find the
alternate allele frequency count for some subset of the samples.  GQT allows
users to identify the target set through a query interface.  This assumes that
the user has a tab-separated file that describes the samples in their BCF/VCF
file, where the each column is a different field, the first row names
the fields, and the remaining rows correspond to the samples in the
BCF/VCF file.  These rows must follow the same order as the samples in the
BCF/VCF file, and there must be the same number of individuals in this file as
there are samples in the BCF/VCF.  This format matches that of the popular PED
format, but GQT does not require any specific fields.  If this file does not already exist, a simple one can be easily generated with the following command:

    # create a header row
    $ echo "Sample Name" > simple.ped
    # extract sample name from the VCF, one per line
    $ bcftools view -h 1kg.chr22.vcf.gz \
        | grep "^#CHROM" \
        | cut -f 10- \
        | tr '\t' '\n' \
        >> simple.ped
    $ head simple.ped
    Sample Name
    HG00096
    HG00097
    HG00099
    HG00100
    HG00101
    HG00102
    HG00103
    HG00104
    HG00106
    $ ~/src/gqt/bin/gqt convert ped -i simple.ped

This convert command will created a database named `simple.ped.db`.  With this
database, which consists of only the sample name (NOTE: spaces in the field
name are converted to underscores("\_")), we can run queries that included
specific individuals

    $ gqt query \
        -i 1kg.chr22.vcf.gz.gqt \
        -d simple.ped.db \
        -p "Sample_Name = 'HG00096' OR Sample_Name = 'HG00097'" \
        -g "count(HET)" \
        > HG00096_7.count_het.vcf

*Step 4*. More sophisticated queries are possible when a more extensive data
base is available.  1000 genomes has a panel file that can be used with some
slight modifications.

    # download the file
    $ wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel
    # add header line, and collapse extra column
    $ (echo -e "Sample ID\tPopulation\tSuper Population\tPlatform"; cat phase1_integrated_calls.20101123.ALL.panel | sed -e "s/,\t/,/") \
        > phase1.panel
    $ head phase1.panel
    Sample ID   Population  Super Population    Platform
    HG00096 GBR EUR ILLUMINA
    HG00097 GBR EUR ABI_SOLID
    HG00099 GBR EUR ABI_SOLID
    HG00100 GBR EUR ILLUMINA
    HG00101 GBR EUR ABI_SOLID,ILLUMINA
    HG00102 GBR EUR ABI_SOLID,ILLUMINA
    HG00103 GBR EUR ILLUMINA
    HG00104 GBR EUR ABI_SOLID
    HG00106 GBR EUR ABI_SOLID,ILLUMINA

*Step 5* Covert the panel file to a db.  This command will create a file called
`phase1.panel.db`.  You can specify the output file name by using the `-o`
option.

    $ gqt convert ped \
        -i phase1.panel

*Step 6* Construct a more interesting query based on the panel file. In this
case, get the count of the non-ref alleles observed for each in the GBR
population and the percent of non-ref alleles in the FIN population.

    $ gqt query \
         -i 1kg.chr22.vcf.gz.gqt \
         -d phase1.panel.db \
         -p "Population=='GBR'" \
         -g "count(HET HOMO_ALT)" \
         -p "Population=='FIN'" \
         -g "pct(HET HOMO_ALT)" \
      > nonref_count_GBR_pct_FIN.vcf


Output
=============

By default, GQT returns valid VCF with results appended to the INFO field.  For
the query:

    $ gqt query \
         -i 1kg.chr22.vcf.gz.gqt \
         -d phase1.panel.db \
         -p "Population=='GBR'" \
         -g "count(HET HOMO_ALT)" \
         -p "Population=='FIN'" \
         -g "pct(HET HOMO_ALT)" \
      > nonref_count_GBR_pct_FIN.vcf

Two lines are added to header:

    ##INFO=<ID=GTQ_0,Number=1,Type=Integer,Description="GQT count result from query 0">
    ##INFO=<ID=GTQ_1,Number=1,Type=Float,Description="GQT percent result from query 1">

Where query 0 is: `-p "Population=='GBR'" -g "count(HET HOMO_ALT)"` and query 1
is `-p "Population=='FIN'" -g "pct(HET HOMO_ALT)"`.

Each variant in `nonref_count_GBR_pct_FIN.vcf` contains two extra values in the
`INFO` field that correspond to the queries.  For example,
`GTQ_0=2;GTQ_1=0.010753` indication that 2 individuals from the GBR population
and 0.01075 percent of the FIN population had a non-ref allele for that variant.

Filtering with `bcftools`
=============

GQT returns valid VCF, which allows for further filtering using `bcftools`.  For
example, to get just the variants that :
- have an average posterior probability of at least 0.99, 
- at least 10 GRB individuals have a non-ref allele
- more than 10% of the FIN population have a non-ref allele

    $ gqt query \
         -i 1kg.chr22.vcf.gz.gqt \
         -d phase1.panel.db \
         -p "Population=='GBR'" \
         -g "count(HET HOMO_ALT)" \
         -p "Population=='FIN'" \
         -g "pct(HET HOMO_ALT)" \
    | bcftools view -i 'AVGPOST>0.99 && GTQ_0>10 && GTQ_1>0.1'
