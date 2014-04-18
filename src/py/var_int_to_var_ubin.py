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

parser.add_option("-b",
    "--ubin_file",
    dest="ubin_file",
    help="Uncompressed variant binary output file")

parser.add_option("-p",
    "--print",
    action="store_true", default=False,
    dest="print_to_screen",
    help="Print ints to screen")


(options, args) = parser.parse_args()

if not options.var_int_file:
    parser.error('Variant int file not given')

if not options.print_to_screen and not options.ubin_file:
    parser.error('Uncompressed varaint binary output file not given')

f = open(options.var_int_file,'r')

if options.ubin_file:
    f_out = open(options.ubin_file, 'wb')

tot_vars = -1

 

if options.ubin_file:
    data = array.array('I')
    data.append(1)


for l in f:
    output_int = 0
    num_vars = 0
    A = l.rstrip().split(' ')
    if tot_vars == -1:
        tot_vars = len(A)
        if options.print_to_screen:
            print tot_vars
        else:
            data.append(tot_vars)
    
    for a in A:
        output_int |= int(a) << (30 - num_vars * 2)
        num_vars += 1
        if num_vars == 16:
            if options.print_to_screen:
                print output_int,
            else:
                data.append(output_int)
            output_int = 0
            num_vars = 0
    if num_vars > 0:
        if options.print_to_screen:
            print output_int,
        else:
            data.append(output_int)
    if options.print_to_screen:
        print
    else:
        data.tofile(f_out)
        data = array.array('I')
f.close()
if not options.print_to_screen:
    f_out.close()

