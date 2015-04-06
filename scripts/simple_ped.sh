#!/bin/bash

############################################################
#  Program:
#  Author :
############################################################


## BEGIN SCRIPT
usage()
{
    cat << EOF

usage: $0 OPTIONS

OPTIONS can be:
    -h      Show this message
    -f      BCF/VCF/compressed VCF file
    -o      output file

EOF
}

# Show usage when there are no arguments.
if test -z "$1"
then
    usage
    exit
fi

VERBOSE=
FILENAME=
OUT=''

# Check options passed in.
while getopts "h f:o:" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        f)
            FILENAME=$OPTARG
            ;;
        o)
            OUT=$OPTARG
            ;;
        ?)
            usage
            exit
            ;;
    esac
done


if [ -z "$OUT" ]
then
    echo "Sample Name"
    bcftools view -h $FILENAME \
        | grep "^#CHROM" \
        | cut -f 10- \
        | tr '\t' '\n' 
else
    echo "Sample Name" > $OUT
    bcftools view -h $FILENAME \
        | grep "^#CHROM" \
        | cut -f 10- \
        | tr '\t' '\n' \
    >> $OUT
fi
