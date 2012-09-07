"""
Work in progress . . . use stuff in here to build proper tests . . .
"""

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib
from scurgen import hilbert
import numpy as np

dim = 16
length = 100000

h = hilbert.HilbertNormalized(dim, length)

cells_to_fill = 5
dist_per_cell = float(length) / (dim * dim)
dist = cells_to_fill * dist_per_cell - 1

# Do the actual work
h.update(0, dist)
h.update(50 * dist_per_cell, 75 * dist_per_cell, value=0.5)



# Plot everybody
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111)
mappable = ax.imshow(
    h.matrix, interpolation='nearest', cmap=matplotlib.cm.Spectral_r)
x, y, labels = h.curve()
x = x.astype(float)
y = y.astype(float)
LIMIT = 120
for xi, yi, li in zip(x, y, labels)[:LIMIT]:
    ax.text(xi, yi, li, size=8)
ax.plot(x, y, color='k')
plt.colorbar(mappable)
ax.axis('tight')

from scipy.spatial import Delaunay

xi = x[50:75].astype(int)
yi = y[50:75].astype(int)

ax.scatter(xi, yi, c=range(len(xi)), s=50, zorder=100)

# Experimenting with convex hulls to eventually plot chromosome boundaries.
# But really what we need is a concave hull...

d = Delaunay(zip(xi, yi))
edges = set()
edge_points = []
def add_edge(i, j):
    if (i, j) in edges or (j, i) in edges:
        return
    edges.add((i, j))
    edge_points.append(d.points[[i, j]])

for ia, ib in d.convex_hull:
    add_edge(ia, ib)
lines = LineCollection(edge_points, linewidth=10)
ax.add_collection(lines)

m = np.zeros_like(h.matrix)

plt.show()

