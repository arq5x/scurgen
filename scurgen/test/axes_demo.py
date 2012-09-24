import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gs
import numpy as np

fig = plt.figure(figsize=(10,10))

nchroms = 21
nrows = int(np.round(np.sqrt(nchroms)))
ncols = nrows

nfiles = 3

CHROM = dict(
    left= 0.05,
    right=0.8,
    top=0.9,
    bottom=0.2,
    wspace=0.1,
    hspace=0.1)

SLIDER_PAD = 0.01
SLIDER = dict(
    left=CHROM['left'],
    right=CHROM['right'],
    bottom=0.1,
    top=CHROM['bottom'] - SLIDER_PAD,
    hspace=0.5,
    )

CBAR_PAD = 0.01
CBAR = dict(
    left=CHROM['right'] + CBAR_PAD,
    right=0.9,
    wspace=SLIDER['hspace'],
    top=CHROM['top'],
    bottom=CHROM['bottom'],
    )

CHECKS = dict(
    top=SLIDER['top'],
    bottom=SLIDER['bottom'],
    left=SLIDER['right'] + CBAR_PAD,
    right=CBAR['right'],
    wspace=CBAR['wspace'],
    hspace=SLIDER['hspace'])

chroms = gs.GridSpec(nrows, ncols)
chroms.update(**CHROM)
axs1 = [plt.subplot(i) for i in chroms]

sliders = gs.GridSpec(nfiles, 1)
sliders.update(**SLIDER)
axs2 = [plt.subplot(i) for i in sliders]

colorbars = gs.GridSpec(1, nfiles)
colorbars.update(**CBAR)
axs3 = [plt.subplot(i) for i in colorbars]

checks = gs.GridSpec(nfiles, nfiles)
checks.update(**CHECKS)
axs4 = [plt.subplot(checks[i, i]) for i in range(nfiles)]
for ax in axs2 + axs4:
    ax.set_xticks([])
    ax.set_yticks([])
plt.show()
