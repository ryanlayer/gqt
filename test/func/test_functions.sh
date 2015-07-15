#!/bin/bash

############################################################
#  Program: test_functions
#  Author : Ryan M Layer ryan.layer@gmail.com
############################################################

BCFTOOLS=bcftools
GQT=../../bin/gqt
DATA_PATH=../data
BCF=$DATA_PATH/10.1e4.var.bcf

TMP_O=$TMPDIR/o
TMP_E=$TMPDIR/e
OUTVAL=
ERRVAL=
RETVAL=
CMD=
VERBOSE=


#{{{ Command line parsing
usage()
{
    cat << EOF

usage: $0 OPTIONS

OPTIONS can be:
    -h      Show this message
    -v      Print success messages
EOF
}

# Check options passed in.
while getopts "h v" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
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
#}}}

#{{{ support fuctions
#{{{ exit codes
EX_OK=0

#The command was used incorrectly, e.g., with the wrong number of arguments, a
#bad flag, a bad syntax in a parameter, or whatever.
EX_USAGE=64

#The input data was incorrect in some way.  This should only be used for user's
#data and not system files.
EX_DATAERR=65

#An input file (not a system file) did not exist or was not readable.  This
#could also include errors like ``No message'' to a mailer (if it cared to
#catch it).
EX_NOINPUT=66

#The user specified did not exist.  This might be used for mail addresses or
#remote logins.
EX_NOUSER=67

#The host specified did not exist.  This is used in mail addresses or network
#requests.
EX_NOHOST=68

#A service is unavailable.  This can occur if a support program or file does
#not exist.  This can also be used as a catchall message when something you
#wanted to do doesn't work, but you don't know why.
EX_UNAVAILABLE=69

#An internal software error has been detected.  This should be limited to
#non-operating system related errors as possible.
EX_SOFTWARE=70

#An operating system error has been detected.  This is intended to be used for
#such things as ``cannot fork'', ``cannot create pipe'', or the like.  It
#includes things like getuid returning a user that does not exist in the passwd
#file.
EX_OSERR=71

#Some system file (e.g., /etc/passwd, /var/run/utmp, etc.) does not exist,
#cannot be opened, or has some sort of error (e.g., syntax error).
EX_OSFILE=72

#A (user specified) output file cannot be created.
EX_CANTCREAT=73

#An error occurred while doing I/O on some file.
EX_IOERR=74

#Temporary failure, indicating something that is not really an error.  In
#sendmail, this means that a mailer (e.g.) could not create a connection, and
#the request should be reattempted later.
EX_TEMPFAIL=75

#The remote system returned something that was ``not possible'' during a
#protocol exchange.
EX_PROTOCOL=76

#You did not have sufficient permission to perform the operation.  This is not
#intended for file system problems, which should use EX_NOINPUT or
#EX_CANTCREAT, but rather for higher level permissions.
EX_NOPERM=77

#Something was found in an unconfigured or misconfigured state.
EX_CONFIG=78
#}}}

#{{{ function run {
function run {
    CMD=$@
    $CMD > $TMP_O 2> $TMP_E
    RETVAL=$?
    OUTVAL=`cat $TMP_O`
    ERRVAL=`cat $TMP_E`
    rm $TMP_O $TMP_E
}
#}}}

#{{{ function print_exit_code {
function print_exit_code {
    case $1 in
        $EX_OK)
            echo "EX_OK"
            ;;
        $EX_USAGE)
            echo "EX_USAGE"
            ;;
        $EX_DATAERR)
            echo "EX_DATAERR"
            ;;
        $EX_NOINPUT)
            echo "EX_NOINPUT"
            ;;
        $EX_NOUSER)
            echo "EX_NOUSER"
            ;;
        $EX_NOHOST)
            echo "EX_NOHOST"
            ;;
        $EX_UNAVAILABLE)
            echo "EX_UNAVAILABLE"
            ;;
        $EX_SOFTWARE)
            echo "EX_SOFTWARE"
            ;;
        $EX_OSERR)
            echo "EX_OSERR"
            ;;
        $EX_OSFILE)
            echo "EX_OSFILE"
            ;;
        $EX_CANTCREAT)
            echo "EX_CANTCREAT"
            ;;
        $EX_IOERR)
            echo "EX_IOERR"
            ;;
        $EX_TEMPFAIL)
            echo "EX_TEMPFAIL"
            ;;
        $EX_PROTOCOL)
            echo "EX_PROTOCOL"
            ;;
        $EX_NOPERM)
            echo "EX_NOPERM"
            ;;
        $EX_CONFIG)
            echo "EX_CONFIG"
            ;;
    esac
}
#}}}

#{{{function assert_exit_code {
function assert_exit_code {
    E=$(print_exit_code $1)
    O=$(print_exit_code $RETVAL)
    if [ $RETVAL -ne $1 ]
    then
        echo -e "FAILURE EXIT CODE($2): \"$CMD\""
        echo -e "-->\texpected $E, observed $O"
        exit
    else
        echo -e "SUCCESS EXIT CODE($2): \"$CMD\""
        if [ $VERBOSE ] 
        then
            echo -e "-->\texpected $E, observed $O"
        fi
    fi
}
#}}}

#{{{ function assert_no_stdout {
function assert_no_stdout {
    if [ -n "$OUTVAL" ]
    then
        echo -e "FAILURE NON-EMPTY STDOUT($1): \"$CMD\""
        echo -e "-->\t$OUTVAL"
        exit
    else
        echo -e "SUCCESS NON-EMPTY STDOUT($1): \"$CMD\""
    fi
}
#}}}

#{{{function assert_stderr {
function assert_stderr {
    if [ -z "$ERRVAL" ]
    then
        echo -e "FAILURE EMPTY STDERR($1): \"$CMD\""
        exit
    else
        echo -e "SUCESS EMPTY STDERR($1): \"$CMD\""
        if [ $VERBOSE ] 
        then
            echo -e "-->\t$ERRVAL"
        fi
    fi
}
#}}}

#{{{ function assert_fail_to_stderr {
function assert_fail_to_stderr {
    assert_exit_code $1 $2
    assert_no_stdout $2
    assert_stderr $2
}
#}}}

#{{{ function make_index {
function make_index {
    $BCFTOOLS index -f $BCF
    $GQT convert bcf \
        -i $BCF \
        2> /dev/null
    $GQT convert ped \
        -i $BCF 2>/dev/null
}
#}}}

#{{{ function rm_index {
function rm_index {
    rm $BCF.vid
    rm $BCF.gqt
    rm $BCF.db
    rm $BCF.csi
    rm $BCF.bim
}
#}}}
#}}}
