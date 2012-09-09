from scurgen.hilbert import d2xy
import numpy as np
from matplotlib import pyplot as plt


def annotated_hilbert(dim, dist_limit=None, ax=None):
    """
    Plots a Hilbert curve that is `dim` cells on a side.

    Labels the first `dist_limit` points.

    If `ax` is provided, plot on that; otherwise plot on a new figure.
    """
    if not ax:
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)
    xs, ys = [], []
    if dist_limit is None:
        dist_limit = dim * dim
    for i in range(dim * dim):
        x, y = d2xy(dim, i)
        ix = y
        iy = dim - x - 1
        xs.append(ix)
        ys.append(iy)
        if i < dist_limit:
            ax.text(ix, iy, '%s' % i, size=10)
    ax.plot(xs, ys, color='0.5', linewidth=0.5)
    ax.scatter(xs, ys, marker='o', c=range(dim * dim), s=10, linewidths=0)
    ax.axis(xmin=-0.5, xmax=dim - 0.05, ymin=-0.5, ymax=dim - 0.5)
    ax.set_title('Hilbert curve, order=%d' % (np.log2(dim)))

if __name__ == "__main__":
    fig = plt.figure(figsize=(17, 10))
    rows = 2
    cols = 3
    i = 1
    for row in range(rows):
        for col in range(cols):
            dim = 2 ** i
            ax = plt.subplot2grid((rows, cols), (row, col))
            annotated_hilbert(dim, dist_limit=10, ax=ax)
            i += 1
    plt.show()
