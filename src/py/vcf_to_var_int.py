#!/usr/bin/env python
import sys
import numpy as np

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-v",
    "--vcf_file",
    dest="vcf_file",
    help="VCF file")

(options, args) = parser.parse_args()

if not options.vcf_file:
    parser.error('VCF file not given')

f = open(options.vcf_file,'r')

for l in f:
    if l[0] != '#':
        A = l.rstrip().split('\t')
        O =[]
        for a in A[9:]:
            g = a.split(':')[0]
            if g == '0|0':
                O.append('0')
            elif g == '1|0' or g == '0|1':
                O.append('1')
            elif g == '1|1':
                O.append('2')
            else:
                O.append('3')

        print " ".join(O)
f.close()
