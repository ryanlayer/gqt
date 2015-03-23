#!/usr/bin/env python
import sys
import numpy as np

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-m",
    "--matrix",
    dest="matrix",
    help="Matrix file")

(options, args) = parser.parse_args()

f = []
if options.matrix:
    f = open(options.matrix,'r')
else:
    f = sys.stdin

M = []
for l in f:
    A = l.rstrip().split('\t')
    M.append(A)

M_l = len(M)
for i in range(M_l):
    o = []
    fill_l = M_l - len(M[i])
    if (M[i][0] == ''):
        fill_l = M_l
    for j in range(fill_l - 1):
        o.append(M[j][i-j-1])
    o.append('0')
    o+=M[i]
    print '\t'.join(o)

f.close()
