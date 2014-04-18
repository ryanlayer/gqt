#!/usr/bin/env python
import sys
import numpy as np
import array

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-v",
    "--var_int_file",
    dest="var_int_file",
    help="Variant int file")

(options, args) = parser.parse_args()

if not options.var_int_file:
    parser.error('Variant int file not given')

f = open(options.var_int_file,'r')

I = []

for l in f:
    output_int = 0
    num_vars = 0
    A = l.rstrip().split()

    if len(I) == 0:
        I = [[] for x in range(len(A))]

    for i in range(len(A)):
        I[i].append(A[i])

for i in I:
    print ' '.join(i)

