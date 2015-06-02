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
rcParams['legend.numpoints'] = 1

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

cor_mat = np.corrcoef(X.T)
eig_val_cor, eig_vec_cor = np.linalg.eig(cor_mat)

eig_pairs_cor = [(np.abs(eig_val_cor[i]), eig_vec_cor[:,i]) for i in range(len(eig_val_cor))]
eig_pairs_cor.sort()
eig_pairs_cor.reverse()
PCA_1 = 0
PCA_2 = 1
print eig_pairs_cor[PCA_1][0]
print eig_pairs_cor[PCA_2][0]

#exit(1)

matrix_w_cor = np.hstack((eig_pairs_cor[PCA_1][1].reshape(len(M),1), eig_pairs_cor[PCA_2][1].reshape(len(M),1)))

X_transf = matrix_w_cor.T.dot(X_std.T).T




#sklearn_pca = PCA(n_components=4)
#X_transf = sklearn_pca.fit_transform(X_std)
#
#if (options.mirror_m):
    #X_transf = X_transf * -1

X_transf[:,0] = X_transf[:,0] * -1
#X_transf[:,1] = X_transf[:,1] * -1

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
#colors = ['cyan', 'lightgreen', 'lightcoral', 'violet', 'grey']
color_map = {}
for i in range(len(U)):
    color_map[U[i]] = colors[i]

print color_map

matplotlib.rcParams.update({'font.size': 12})
#fig = plt.figure(figsize=(5,5), dpi=300)
#fig = matplotlib.pyplot.figure(figsize=(5,5),dpi=300,facecolor='black')
fig = matplotlib.pyplot.figure(figsize=(5,5),dpi=300)
ax = fig.add_subplot(1,1,1)
#ax = fig.add_subplot(1,1,1,axisbg='k')
ax.tick_params(labelsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#ax.spines['bottom'].set_color('white')
#ax.spines['left'].set_color('white')
#ax.title.set_color('white')
#ax.yaxis.label.set_color('white')
#ax.xaxis.label.set_color('white')
#ax.tick_params(axis='x', colors='white')
#ax.tick_params(axis='y', colors='white')
#ax.set_xlabel("Principal component 1")
#ax.set_ylabel("Principal component 2")
ax.set_xlabel("EV = " + str(round(eig_pairs_cor[PCA_1][0],2)))
ax.set_ylabel("EV = " + str(round(eig_pairs_cor[PCA_2][0],2)))

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

#plt.legend(frameon=False, fontsize=10,labelspacing=0.25,numpoints=1)
#plt.savefig(options.output_file,bbox_inches='tight')
l1=ax.legend(loc='lower right', \
             labelspacing=0.25,\
             frameon=False, \
             fontsize=12, \
             ncol=1)

yticks, yticklabels = matplotlib.pyplot.yticks()
ymin = (3*yticks[0] - yticks[1])/2.
ymax = (3*yticks[-1] - yticks[-2])/2.
matplotlib.pyplot.ylim(ymin, ymax)
matplotlib.pyplot.yticks(yticks)

xticks, xticklabels = matplotlib.pyplot.xticks()
xmin = (3*xticks[0] - xticks[1])/2.
xmax = (3*xticks[-1] - xticks[-2])/2.
matplotlib.pyplot.xlim(xmin, xmax)
matplotlib.pyplot.xticks(xticks)



#for text in l1.get_texts():
    #matplotlib.pyplot.setp(text,color='white')

####################################################
matplotlib.pyplot.savefig(options.output_file,\
                          bbox_inches='tight')

