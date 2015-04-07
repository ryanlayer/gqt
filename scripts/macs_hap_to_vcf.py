#!/usr/bin/env python
import sys

print '##fileformat=VCFv4.1"'
print '##FILTER=<ID=PASS,Description="All filters passed">'
print '##INFO=<ID=N,Type=String>'
print '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
print '##contig=<ID=1,length=3000000000>'
print '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + \
        '\t'.join(["I" + str(x) for x in range(int(sys.argv[1]))])
pos=1
for l in sys.stdin:
    if l[0:4] == 'SITE':
        A = l.rstrip().split()
        H = A[4]
        O = []
        for i in range(0,len(H),2):
            #if len(O) != 0:
                #O += ' '
            #O += '\t' + H[i] + '|' + H[i+1]
            O.append(H[i] + '|' + H[i+1])
        print '\t'.join(['1', \
                         str(pos), \
                         'V' + str(pos), \
                         'A', \
                         'T', \
                         '100', \
                         'PASS', \
                         'N=A', \
                         'GT'] + O)
        pos += 1
