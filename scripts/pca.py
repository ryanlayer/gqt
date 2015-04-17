#!/usr/bin/env python
import sys
from operator import add
import numpy as np
from sklearn.decomposition import PCA
from sklearn import preprocessing
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-c",
    dest="cat_file",
    help="Categories")

parser.add_option("-o",
    dest="output_file",
    help="Output file")


(options, args) = parser.parse_args()
if not options.output_file:
    parser.error('Output file not given')

M = []

for l in sys.stdin:
    A = [float(x) for x in l.rstrip().split()]
    M.append(A)

X = np.array(M)
X_std = preprocessing.StandardScaler().fit_transform(X)
sklearn_pca = PCA(n_components=2)
X_transf = sklearn_pca.fit_transform(X_std)


#print X_transf[:,0]
#print X_transf[:,1]

f = open(options.cat_file)
C=[]
for l in f:
    C.append(l.rstrip())
f.close()

tmp=set(C)
U=[]
for t in tmp:
    U.append(t)
    
colors = cm.rainbow(np.linspace(0, 1, len(U)))
color_map = {}
for i in range(len(U)):
    color_map[U[i]] = colors[i]

print color_map

fig = plt.figure(figsize=(5,5), dpi=300)
ax = fig.add_subplot(1,1,1)

for cat in U:
    print cat
    for idx in [i for i, x in enumerate(C) if x == cat]:
        plt.scatter(X_transf[idx,0], \
                    X_transf[idx,1], \
                    c=color_map[cat],
                    label=cat)

#plt.legend()
plt.savefig(options.output_file,bbox_inches='tight')

