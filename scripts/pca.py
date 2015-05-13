#!/usr/bin/env python
import sys
from operator import add
import numpy as np
from sklearn.decomposition import PCA
from sklearn import preprocessing
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-c",
    dest="cat_file",
    help="Categories")

parser.add_option("-o",
    dest="output_file",
    help="Output file")

parser.add_option("-m",
                  "--mirror",
                  action="store_true",
                  dest="mirror_m",
                  default=False)
                  
(options, args) = parser.parse_args()
if not options.output_file:
    parser.error('Output file not given')

if not options.cat_file:
    parser.error('Category file not given')

M = []

for l in sys.stdin:
    A = [float(x) for x in l.rstrip().split()]
    M.append(A)

X = np.array(M)
X_std = preprocessing.StandardScaler().fit_transform(X)
sklearn_pca = PCA(n_components=4)
X_transf = sklearn_pca.fit_transform(X_std)

if (options.mirror_m):
    X_transf = X_transf * -1

X_transf[:,0] = X_transf[:,0] * -1

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

matplotlib.rcParams.update({'font.size': 12})
fig = plt.figure(figsize=(5,5), dpi=300)
ax = fig.add_subplot(1,1,1)
ax.tick_params(labelsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
l = []
for cat in U:
    print cat
#    tmp_X = []
#    tmp_Y = []
#    for idx in [i for i, x in enumerate(C) if x == cat]:
#        tmp_X.append(X_transf[idx,0])
#        tmp_Y.append(X_transf[idx,1])
#    print tmp_X
#    print tmp_Y
#    plt.scatter(tmp_X, tmp_Y, X_transf[idx,1], c=color_map[cat], label=cat)

    idxs = [i for i, x in enumerate(C) if x == cat]

    ax.plot(X_transf[idxs[0],0], \
                X_transf[idxs[0],1], \
                'o',
                c=color_map[cat],
                label=cat)

    for idx in idxs[1:]:
        ax.plot(X_transf[idx,0], \
                X_transf[idx,1], \
                'o',
                c=color_map[cat])

plt.legend(frameon=False, fontsize=10,labelspacing=0.25,numpoints=1)
plt.savefig(options.output_file,bbox_inches='tight')

