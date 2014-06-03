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
    -i      input file
    -o      output file 
    -v      Verbose (boolean)

EOF
}

# Show usage when there are no arguments.
if test -z "$1"
then
    usage
    exit
fi

VERBOSE=
INFILE=
OUTFILE=

# Check options passed in.
while getopts "h i:o:v" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        i)
            INFILE=$OPTARG
            ;;
        o)
            OUTFILE=$OPTARG
            ;;
        v)
            VERBOSE=1
            ;;
        ?)
            usage
            exit
            ;;
    esac
done

# Do something with the arguments...

## END SCRIPT

FIELDS=`head -n 1 $INFILE | awk '{print NF;}'`
RECORDS=`cat $INFILE | wc -l`

echo $FIELDS > $OUTFILE
echo $RECORDS >> $OUTFILE
cat $INFILE >> $OUTFILE
