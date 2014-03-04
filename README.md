Encoded files are defined by their rows: variants (var) or individuals (ind)
and by the encoding scheme: integer (i), uncompressed binary (ub)

Variant integer representation (vari):

This is space-seperated, ASCI text file where each line gives the integer
representation for the genotypes for a set of individuals where

Homozygous Ref: 0
Heterozygous:   1
Homozygous Alt: 2
Other:          3


Variant uncompressed binary interger representation (varub):

This file is writen in binary as a series of 32-bit ints.  The first int is 1
and verifies the file type, the second int gives the number of individuals, and
the following ints (some multiple of the number of individuals) encode the
genotype of some number of variants.  Each int encodes at most 16 genotypes
using 2 bits per genotype where:

Homozygous Ref: 00
Heterozygous:   01
Homozygous Alt: 10
Other:          11

When the number of individuals is not a multiple of 16, the last int is padded
with zeros.
