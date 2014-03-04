#!/usr/bin/env python
import sys
import struct
from optparse import OptionParser


def var_ubin_to_var_int(var_ubin):
    V = []
    for i in range(30,-1,-2):
        v = (var_ubin >> i) & 3
        V.append(v)
    return V

parser = OptionParser()

parser.add_option("-b",
    "--ubin_file",
    dest="ubin_file",
    help="Uncompressed variant binary file")

(options, args) = parser.parse_args()

if not options.ubin_file:
    parser.error('Binary file not given')

f = open(options.ubin_file,'rb')

vars_per_int = 16

validate = -1
num_inds = -1
num_ints = -1
col = []

try:
    byte = f.read(4)

    while byte != "":

        curr_int = struct.unpack('I', byte)[0]

        if validate == -1:
            if curr_int == 1:
                validate = 1
            else:
                exit(1)
        elif num_inds == -1:
            num_inds = curr_int
            num_ints = int((num_inds + vars_per_int - 1)/vars_per_int)
        else:
            col.append(curr_int)

            if len(col) == num_ints:
                V = []
                for c in col:
                    V += var_ubin_to_var_int(c)
                R = V[:num_inds]
                print " ".join([str(r) for r in R])
                col = []

        byte = f.read(4)
finally:
    f.close()
