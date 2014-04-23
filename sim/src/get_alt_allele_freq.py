#!/usr/bin/env python
import sys

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-d",
    "--data_file",
    dest="data_file",
    help="Data file")

(options, args) = parser.parse_args()

if not options.data_file:
    parser.error('Data file not given')

f = open(options.data_file,'r')


for l in f:
    tot = 0
    A = l.rstrip().split()

    for a in A:
        if int(a) > 0:
            tot += 1
    print tot

f.close()
