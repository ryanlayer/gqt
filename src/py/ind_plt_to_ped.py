#!/usr/bin/env python
import sys
import numpy as np

from optparse import OptionParser

parser = OptionParser()

parser.add_option('-p',
    '--plt',
    dest='plt_file',
    help='plain text file')

(options, args) = parser.parse_args()

if not options.plt_file:
    parser.error('Plain text file not given')

f = open(options.plt_file,'r')

skip = int(f.readline())
skip = int(f.readline())


j = 1
for l in f:
    O = []

    for i in l.rstrip().split():
        if i == '0':
            O.append('1')
            O.append('1')
        elif i == '1':
            O.append('1')
            O.append('2')
        elif i == '2':
            O.append('2')
            O.append('2')
        elif i == '3':
            O.append('0')
            O.append('0')

    H = ['I'+str(j), '1', '0', '0' , '1' , '1']
    
    print ' '.join(H+O)

    j+=1
