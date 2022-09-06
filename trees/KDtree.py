
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree

pts = [[1,6,7,2],
       [2,7,8,3],
       [3,8,9,4],
       [4,9,10,5],
       [9,12,11,10],
       [8,13,12,9],
       [7,14,13,8],
       [6,15,14,7],
       [15,16,17,14],
       [14,17,18,13],
       [13,18,19,12],
       [12,19,20,11]]

kdtree = KDTree(pts)
prs = kdtree.query_pairs(r=2)

print("prs: ", prs)

#plt.figure(figsize=(12, 4))
#plt.plot(pts[:, 0], pts[:, 1], "xk", markersize=14)

#for (i, j) in prs:
#    plt.plot([pts[i, 0], pts[j, 0]],
#            [pts[i, 1], pts[j, 1]], "-r")
#plt.show()
