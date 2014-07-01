#!/usr/bin/env python
import sys
import numpy as np

from optparse import OptionParser

parser = OptionParser()

parser.add_option('-p',
    '--plt',
    dest='plt_file',
    help='plain text file')

parser.add_option('-v',
    '--var',
    dest='var',
    action="store_true",
    default=False,
    help='by-variant')

parser.add_option('-i',
    '--ind',
    dest='ind',
    action="store_true",
    default=False,
    help='by-individual')

(options, args) = parser.parse_args()

if not options.plt_file:
    parser.error('Plain text file not given')

if (not options.ind) and (not options.var):
    parser.error('Must set either by-variant or by-individual')

if options.ind and options.var:
    parser.error('Must set either by-variant or by-individual, not both')

f = open(options.plt_file,'r')

num_var = int(f.readline())


if options.var:
    num_var = int(f.readline())

#The cols are: chrom, id, genetic distance, and bas-pair position
for i in range(num_var):
    print 1, "v" + str(i), 0, i+1
f.close()
