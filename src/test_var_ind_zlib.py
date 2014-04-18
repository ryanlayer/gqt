#!/usr/bin/env python
import sys
import array
import zlib
import numpy as np

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

var_sizes = []

for l in f:
    no_space = l.rstrip().translate(None, ' ')
    var_sizes.append(len(zlib.compress(no_space,9)))

print len(var_sizes)
print 'median',np.median(var_sizes)
print 'mean',np.mean(var_sizes)
print 'std',np.std(var_sizes)
print '[',min(var_sizes),',',max(var_sizes),']'
H = np.histogram(var_sizes)
print H[0]
print H[1]
