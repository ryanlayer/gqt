#!/usr/bin/env python
import sys
import genotq
import struct
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
    dest="variant",
    type="int",
    help="Variant to report")

(options, args) = parser.parse_args()

if not options.ubin_file:
    parser.error('Binary file not given')

if not options.variant:
    parser.error('Variant not given')

f = open(options.ubin_file,'rb')

vars_per_int = 16

validate = -1
num_inds = -1
num_ints = -1
col = []

try:
    #validate
    byte = f.read(4)
    curr_int = struct.unpack('I', byte)[0]

    if curr_int != 1:
        exit(1)

    #get size
    byte = f.read(4)
    curr_int = struct.unpack('I', byte)[0]

    num_inds = curr_int
    num_ints = int((num_inds + vars_per_int - 1)/vars_per_int)


    f.seek(8 + ((options.variant - 1) * num_ints)*4 )

    byte = f.read(4)

    while byte != "":
        curr_int = struct.unpack('I', byte)[0]

        col.append(curr_int)

        if len(col) == num_ints:
            V = []
            for c in col:
                V += genotq.var_ubin_to_var_int(c)
            R = V[:num_inds]
            print " ".join([str(r) for r in R])
            break

        byte = f.read(4)
finally:
    f.close()
