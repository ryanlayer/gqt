#!/usr/bin/env python
import sys
import genotq
from optparse import OptionParser


#def var_ubin_to_var_int(var_ubin):
    #V = []
    #for i in range(30,-1,-2):
        #v = (var_ubin >> i) & 3
        #V.append(v)
    #return V

parser = OptionParser()

parser.add_option("-b",
    "--ubin_file",
    dest="ubin_file",
    help="Uncompressed variant binary file")

parser.add_option("-v",
    "--variants",
    dest="variants",
    help="Comma-sepperated list of variants to report")

(options, args) = parser.parse_args()

if not options.ubin_file:
    parser.error('Binary file not given')

if not options.variants:
    parser.error('Variants not given')

[validate,num_inds,num_ints,f] = genotq.var_ubin_validate(options.ubin_file)

if validate != 1:
    print "Incorrect file type"
    exit(1)

R = [int(x) for x in options.variants.split(',')]
V = genotq.var_ubin_get_rows(f, num_ints, R)

for v in V:
    O = []
    for i in v:
        O += genotq.var_ubin_to_var_int(i)

    print " ".join([str(o) for o in O[:num_inds]])

genotq.var_ubin_finalize(f)
