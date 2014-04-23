#!/usr/bin/env python
import sys

for l in sys.stdin:
    if l[0:4] == 'SITE':
        A = l.rstrip().split()
        H = A[4]
        O = ''
        for i in range(0,len(H),2):
            if len(O) != 0:
                O += ' '
            O += str(int(H[i]) + int(H[i+1]))
        print O
