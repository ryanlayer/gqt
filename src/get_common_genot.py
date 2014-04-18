#!/usr/bin/env python
import sys
import numpy as np
import array

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-q",
    "--quiet",
    dest="quiet",
    action="store_true",
    default=False,
    help="Supress output")

parser.add_option("-v",
    "--var_int",
    dest="var_int",
    action="store_true",
    default=False,
    help="Set input to variant int")

parser.add_option("-i",
    "--ind_int",
    dest="ind_int",
    action="store_true",
    default=False,
    help="Set input to individual int")

parser.add_option("-I",
    "--inds",
    dest="inds",
    help="CSV set of individuals to test")

(options, args) = parser.parse_args()

if not options.inds:
    parser.error('Set of individuals to test not given')

if not options.var_int and not options.ind_int:
    parser.error('Must set either variant int or individual int')

if options.var_int and options.ind_int:
    parser.error('Must set either variant int or individual int, not both')

I = sorted([int(x) for x in options.inds.split(',')])

common_var = []



if options.var_int:
    for l in sys.stdin:
        A = l.rstrip().split() 
        r = 3 
        for g in [int(A[i]) for i  in I]:
            r = r & (g>0)
        common_var.append(r)

elif options.ind_int:
    line_num = 0
    curr_I = 0
    for l in sys.stdin:
        if line_num == I[curr_I]:
            A = l.rstrip().split() 
            if len(common_var) == 0:
                common_var = [3] *  len(A)

            for i in range(len(A)):
                common_var[i] = common_var[i] & (int(A[i])>0)

            curr_I += 1
        
        if curr_I == len(I):
            break

        line_num += 1
        
if not options.quiet:
    for i in range(len(common_var)):
        if common_var[i] > 0:
            print i

