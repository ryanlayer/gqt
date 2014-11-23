Overview
========

NOTE: GQT is in an alpha state. 

Genome Query Tools (GQT) is a tool and C API for storing and querying
large-scale genotype data sets like those produced by 1000 Genomes. Genotypes
are represented by compressed bitmap indices, which reduce the storage and
compute burden by orders of magnitude. This index can significantly expand the
capabilities of population-scale analyses by providing interactive-speed
queries to data sets with millions of individuals. An example work flow:

Use GQT to create a sample-based index of a BCF file containing 1092
individuals and 494328 variants.

    gqt convert bcf \
        -r 494328 \
        -f 1092 \
        -i 1kg.chr22.bcf 

The same command works for VCF and compressed VCF

    gqt convert bcf \
        -r 494328 \
        -f 1092 \
        -i 1kg.chr22.vcf 

    gqt convert bcf \
        -r 494328 \
        -f 1092 \
        -i 1kg.chr22.vcf.gz

Create a database of the PED file describing the phenotypes and relationships
of the 1092 samples.

    gqt convert ped \
        -i 1kg.ped 

Find all variants where at least 10 GBR individuals are heterozygous.

    gqt query \
        -i 1kg.chr22.bcf.gqt \
        -d 1kg.ped.db \
        -p "Population = 'GBR'" \
        -g "count(HET) >= 10"

Find all variants where all affected individuals are homozygous.

    gqt query \
        -i 1kg.chr22.bcf.gqt \
        -d 1kg.ped.db \
        -p "Phenotype = 2" \
        -g "HOMO_REF”

Further filter variants by total depth with bcftools.

    gqt query \
        -i 1kg.chr22.bcf.gqt \
        -d 1kg.ped.db \
        -p "Population = 'GBR'" \
        -g "count(HET) >= 10" \
    | bcftools view - -i 'DP>50000'


GQT takes a BCF as input and produces three files, a (very small) compressed
index (.gqt), a compressed summary of the variant metadata (.bim), and a
variant ID file (.vid). This process rotates the data, sorts it by alternate
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
command line argument "-r", or by querying an associated PED file
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


Example usage
=============

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
alternate allele frequency count for some subset of the samples.  The first
command simply makes a list of 100 comma-separated-values between 0 and 99.  In
a real query you would clearly choose the IDs associated with your set of
interest.  We will soon move most of these options to a config file that can be
specified on the command line, which will make selecting large groups much less
error prone.
 
    $ Q=`seq 0 99|tr '\n' ',' | sed -e "s/,$//"`
    $ gqt query \
         -i 1kg.chr22.vcf.gz.gqt \
         -n 100 \
         -r $Q

If you have compiled in AVX2 support, then GQT will automatically use those functions to give you much better performance.

*Step 4*. Download the Phase 1 1000 genomes PED file. NOTE: this is from Phase
3, so it does not match exactly but is useful for the example.  Again, shorten the file name.

    $ wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp//release/20130502/integrated_call_samples.20130502.ALL.ped
    $ mv integrated_call_samples.20130502.ALL.ped 1kg.ped

*Step 5* Covert the PED file to a PED db.  This command will create a file called `1kg.ped.db`.  You can specify the output file name by using the `-o` option.

    $ gqt convert ped-db \
        -i 1kg.ped

*Step 6* Submit the same query as above, but this time using a query based on the PED file.

    $ gqt query \
         -i 1kg.chr22.vcf.gz.gqt \
         -d 1kg.ped.db \
         -p "Ind_ID < 100"
         -g "count(HET)"

Or, construct a more interesting query based on the PED file. In this case,
count the non-ref alleles observed for each variant among CEU females.

    $ gqt query \
         -i 1kg.chr22.vcf.gz.gqt \
         -d 1kg.ped.db \
         -p "Population=='CEU' and Gender==2"
         -q "count(HET HOMO_ALT)"
