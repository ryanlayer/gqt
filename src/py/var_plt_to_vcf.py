#!/usr/bin/env python
import sys
import numpy as np

from optparse import OptionParser

parser = OptionParser()

parser.add_option('-v',
    '--var_plt',
    dest='var_plt_file',
    help='by-variant plain text file')

(options, args) = parser.parse_args()

if not options.var_plt_file:
    parser.error('By-variant plain text file not given')

f = open(options.var_plt_file,'r')

#print the required fields

print '##fileformat=VCFv4.1'
print '##INFO=<ID=N,Type=String>'
print '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
print '##contig=<ID=1,length=3000000000>'



num_records = -1
num_fields = -1
header = 0
pos = 1
for l in f:
    if num_fields == -1:
        num_fields = int(l)
    elif num_records == -1:
        num_records = int(l)
    elif header == 0:
        print '#CHROM\t' + \
              'POS\t' + \
              'ID\t' + \
              'REF\t' + \
              'ALT\t' + \
              'QUAL\t' + \
              'FILTER\t' + \
              'INFO\t' + \
              'FORMAT\t' + \
              '\t'.join([ 'I' + str(i) for i in range(num_fields)])
        header = 1
    else:
        A = l.rstrip().split()

        G = []
        for a in A:
            if a == '0':
                G.append('0|0')
            elif a == '1':
                G.append('0|1')
            elif a == '2':
                G.append('1|1')

        O = ['1', \
             str(pos), \
             'V' + str(pos), \
             'A',
             'T',
             '100',
             'PASS',
             'N=A',
             'GT' ] + G
        
        print '\t'.join(O)
        pos += 1
