#!/usr/bin/env python
import sys
import struct
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-b",
    "--binary_file",
    dest="binary_file",
    help="Binary file")

parser.add_option("-c",
    "--count",
    dest="count",
    type="int",
    default=-1,
    help="Number of lines to read")



(options, args) = parser.parse_args()

if not options.binary_file:
    parser.error('Binary file not given')

f = open(options.binary_file,'rb')
count = options.count 

try:
    byte = f.read(4)
    while byte != "" and count != 0:
        print struct.unpack('I', byte)[0]
        count =- 1
        byte = f.read(4)
finally:
    f.close()
