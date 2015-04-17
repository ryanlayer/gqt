#!/usr/bin/env python
import sys
from operator import add
import numpy as np
from sklearn.decomposition import PCA
from sklearn import preprocessing
import matplotlib.pyplot as plt

M = []

for l in sys.stdin:
    A = [float(x) for x in l.rstrip().split()]
    M.append(A)

X = np.array(M)
X_std = preprocessing.StandardScaler().fit_transform(X)
sklearn_pca = PCA(n_components=2)
X_transf = sklearn_pca.fit_transform(X_std)


print X_transf[:,0]
print X_transf[:,1]

plt.scatter(X_transf[:,0], X_transf[:,1])
plt.show()

