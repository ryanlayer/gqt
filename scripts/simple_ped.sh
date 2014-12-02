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

# Check options passed in.
while getopts "h f:" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        f)
            FILENAME=$OPTARG
            ;;
        ?)
            usage
            exit
            ;;
    esac
done

echo "Sample Name"

bcftools view -h $FILENAME \
    | grep "^#CHROM" \
    | cut -f 10- \
    | tr '\t' '\n' 
