Genome Query Tools (GQT) is a tool and C API for storing and querying
large-scale genotype data sets like those produced by 1000 Genomes.  Genotypes
are represented by compressed bitmap indices, which reduce the storage and
compute burden by orders of magnitude. This index can significantly expand the
capabilities of population-scale analyses by providing interactive-speed
queries to data sets with millions of individuals.

GQT takes a BCF as input and produces two files, a compressed index (.wahbm)
and a summary of the variant metadata (.bim). This process rotates the data,
sorts it by alternate allele frequency, converts it to a bitmap index, then
finally compresses the data using the Word Aligned Hybrid (WAH) encoding
(http://dl.acm.org/citation.cfm?id=1132864).  Rotating the typical
variant-row/individual-column form of a BCF such that individuals become rows
allows for more efficient memory access patterns for individual based queries
(e.g., find the alternate allele frequency count for some set of individuals
within the cohort).  Sorting columns (variants) improves compression.
Converting to a bitmap index allows queries to be resolved using bit-wise
logical operations, which can compare 32 genotypes in a single fast operation.
WAH encoding gives near-optimal compression while allowing bit-wise logical
operations without inflation.

The combine result these steps creates an index that is only a fraction of the
source BCF size, and queries against that index complete in seconds.  In the
case of chr22 from 1000 Genomes, the BCF file is 11G and the WAH index is 42M,
and the alternate allele frequency count for 100 of those individuals can be
found in 0.188 seconds.

Installation
============
GQT depends on htslib and hdf5.

*Step 1*. Install htslib.

    git clone https://github.com/samtools/htslib.git
    cd htslib
    make


*Step 2*. Install hdf5. We recommend downloading one of the statically-linked binary distributions from
here: http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.7/obtain5187.html.

*Step 3*. Modify the GQT Makefile by change the `HTS_ROOT` and `HDF_ROOT` variables in src/c/Makfile to
reflect ther locations.

*Step 4*. Compile GQT

    cd gqt/
    make

*Step 5*. Test GQT

    cd gqt/src/test/unit
    make
    cd ../func
    bash functional_tests.sh


Example usage
=============

In the example below, we demonstrate how to create a word-aligned hybrid index
of a VCF file from the 1000 Genomes Project. We then demonstrate how to use
this index to quickly query variants based on sample genotypes (the query
interface is very limited, but we are actively working on it.). Please note
that the conversion of BCF to a WAH index will take roughly 4 hours, as we have
not yet had time to focus on optimizing or parallelizing this conversion.  If
you are impatient, we have posted all of the converted files described above on
our lab [website](http://quinlanlab.cs.virginia.edu/gqt-example/).

*Step 1*. Download the Phase 1 1000 genomes chr22 file.

	$ wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz

*Step 2*. Use the new (and very nice) version of [bcftools](http://samtools.github.io/bcftools/) to convert the file VCF to BCF.

	$ bcftools view ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz > 1kg.chr22.bcf

*Step 3*. Use GQT to make a word-aligned hybrid index of the genotypes in the BCF file. Currently, one must tell GQT how many variants (`-r`) and how many samples (`-f`) are in the file. Future versions will detect this automatically.  To extract this info, one can use the [bcftools](http://samtools.github.io/bcftools/) `stats` command. 

	# use bcftools stats to get the number of variants and individuals in the BCF file.
	$ bcftools stats 1kg.chr22.bcf | \
	    grep SN | \
	    head -4
	# SN, Summary numbers:
	# SN	[2]id	[3]key	[4]value
	SN	0	number of samples:	1092
	SN	0	number of records:	494328

As you can see, there are 494328 variants (records) and 1092 samples. We
provide this information to the GQT `convert` tool and ask it to make a
word-aligned hybrid index (`bcf-wahbm`) of the BCF file. In adddition, the tool
will create a `.bim` file that stores very basic information about the variants
in the original BCF. The current implementation of GQT reports information from
this smaller `.bim` file for speed and simplicity, but future versions will
reach back into the BCF file and report the full BCF record for each variant
that meets the genotype query criteria given to GQT by the user.

	$ gqt convert bcf-wahbm   \
	    -r 494328             \
	    -f 1092               \
	    -i 1kg.chr22.bcf      \
	    -b 1kg.chr22.bcf.bim  \
	    -o 1kg.chr22.bcf.wahbm

*Step 4*.  Query the BCF file using the WAH index created by GQT.  In this case
we are going to simply find the alternate allele frequency count for some
subset of the samples.  The first command simply makes a list of 100
comma-separated-values between 0 and 99.  In a real query you would clearly
choose the ids associated with your set of interest.  We will soon move most of
these options to a config file that can be specified on the command line, which
will make selecting large groups much less error prone.
 
        $ Q=`seq 0 99|tr '\n' ',' | sed -e "s/,$//"`
        $ gqt sum ipwahbm \
            -b 1kg.chr22.bcf.bim \
            -i 1kg.chr22.bcf.wahbm \
            -n 100 \
            -r $Q

If you have compiled in AVX2 support (uncomment line 7 of the Makefile in gqt/src/c) you can get much better performance by using the "-a" option.

        $ gqt sum ipwahbm \
            -a \
            -b 1kg.chr22.bcf.bim \
            -i 1kg.chr22.bcf.wahbm \
            -n 100 \
            -r $Q
