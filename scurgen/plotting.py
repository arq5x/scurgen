import os
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.widgets import Cursor, Slider, RadioButtons
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scurgen.hilbert import HilbertMatrix


def debug_plot(h):
    """
    Quick plot of a HilbertBase subclass that also labels the first 10 cells
    """
    imshow_kwargs = dict(
        interpolation='nearest',
        cmap=matplotlib.cm.Spectral_r)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    mappable = ax.imshow(h.matrix, **imshow_kwargs)
    plt.colorbar(mappable)
    x, y, labels = h.curve()
    ax.plot(x, y, 'k')
    for i in range(10):
        ax.text(x[i], y[i], labels[i])
    ax.axis('tight')
    plt.show()


class HilbertGUI(object):
    def __init__(self, intervals1, intervals2, **kwargs):
        """
        :param intervals1:
            A file supported by pybedtools (BED, VCF, GTF, etc)

        :param intervals2:
            Another file supporte by pybedtools

        :param kwargs:
            Additional keyword arguments are passed to HilbertMatrix (e.g.,
            m_dim, genome, chrom)
        """
        # TODO: lots of possible configuration here.  Possibly use a YAML
        # config file strategy?

        self.h1 = HilbertMatrix(intervals1, **kwargs)
        self.h2 = HilbertMatrix(intervals2, **kwargs)

        # masked values become 0% alpha all the time
        self.h1.mask_low_values()
        self.h2.mask_low_values()

    def plot(self):
        """
        Does most of the work to set up the figure, axes, and widgets.
        """
        self.fig = plt.figure()
        self.ax = plt.Axes(self.fig, (0.1, 0.15, 0.8, 0.8))

        # from axes_grid toolkit
        # TODO: eventually support n-way comparisons
        divider = make_axes_locatable(self.ax)
        self.cax1 = divider.append_axes('right', size=0.2, pad=0.3)
        self.cax2 = divider.append_axes('right', size=0.2, pad=0.3)

        # Alpha sliders.
        # TODO: eventually n-way comparisons; these should be appended on the
        # bottom of self.ax?
        self.slider_ax1 = plt.Axes(self.fig, (0.3, 0.02, 0.3, 0.02))
        self.slider_ax2 = plt.Axes(self.fig, (0.3, 0.07, 0.3, 0.02))

        # Radio buttons axes
        self.radio_ax = plt.Axes(self.fig, (0.75, 0.02, 0.1, 0.1))

        # Add everyone to the fig
        self.fig.add_axes(self.ax)
        self.fig.add_axes(self.slider_ax1)
        self.fig.add_axes(self.slider_ax2)
        self.fig.add_axes(self.radio_ax)

        # plot the matrices on top of each other
        self.mappable1 = self.ax.imshow(
            self.h1.matrix, interpolation='nearest',
            cmap=matplotlib.cm.Reds)

        self.mappable2 = self.ax.imshow(
            self.h2.matrix, interpolation='nearest', cmap=matplotlib.cm.Blues)

        # Initialize alphas
        self.mappable1.set_alpha(0.5)
        self.mappable2.set_alpha(0.5)

        # colorbars
        self.cbar1 = plt.colorbar(self.mappable1, cax=self.cax1)
        self.cbar2 = plt.colorbar(self.mappable2, cax=self.cax2)

        # formatting: remove ticklabels on the main axes since they're not
        # meaningful
        self.ax.set_xticks([])
        self.ax.set_yticks([])

        # Tweak colorbar labels
        for cbar in [self.cbar1, self.cbar2]:
            for txt in cbar.ax.get_yticklabels():
                txt.set_size(8)

        # Set up sliders with sensible default labels
        self.slider2 = Slider(
            self.slider_ax1,
            os.path.basename(self.h1.file),
            valmin=0,
            valmax=1,
            valinit=0.5)

        self.slider1 = Slider(
            self.slider_ax2,
            os.path.basename(self.h2.file),
            valmin=0,
            valmax=1,
            valinit=0.5)

        # For controlling log/linear color scale
        self.radio = RadioButtons(
            self.radio_ax,
            labels=['log', 'linear'],
            active=1)

        # These colors look nice with the current colormaps (Reds and Blues)
        # but TODO: should eventually be configurable
        self.slider1.poly.set_color('#7a0510')
        self.slider2.poly.set_color('#08316d')

        # Label tweaks
        self.slider1.label.set_size(10)
        self.slider2.label.set_size(10)

        # Callback changes alpha
        self.slider1.on_changed(self._slider_callback_factory(self.mappable1, self.cbar1))
        self.slider2.on_changed(self._slider_callback_factory(self.mappable2, self.cbar2))

        # Callback changes color scale
        self.radio.on_clicked(self._radio_callback)

    def _radio_callback(self, label):
        # Hello? Am I the 9th caller!?
        if label == 'log':
            self._log()
        elif label == 'linear':
            self._linear()
        else:
            raise ValueError("unspecified label for radio button")

    def _slider_callback_factory(self, mappable, cbar):
        """
        Given a mappable (i.e., object returned from imshow()), return
        a function that will modify it based on the slider's value.
        """
        def _slider_callback(x):
            mappable.set_alpha(x)
            cbar.set_alpha(x)
            cbar.update_normal(mappable)

        return _slider_callback

    def _log(self):
        """
        Update colomaps of the plotted images to use log-scaled color
        """
        norm1 = matplotlib.colors.LogNorm(
            vmin=self.h1.matrix.min(), vmax=self.h1.matrix.max())
        norm2 = matplotlib.colors.LogNorm(
            vmin=self.h2.matrix.min(), vmax=self.h2.matrix.max())
        self.mappable1.set_norm(norm1)
        self.mappable2.set_norm(norm2)
        self.cbar1.set_norm(norm1)
        self.cbar2.set_norm(norm2)
        self.cbar1.update_normal(self.mappable1)
        self.cbar2.update_normal(self.mappable2)
        plt.draw()

    def _linear(self):
        """
        Update colormaps of the plotted images to use linear-scaled color
        """
        norm1 = matplotlib.colors.Normalize(
            vmin=self.h1.matrix.min(), vmax=self.h1.matrix.max())
        norm2 = matplotlib.colors.Normalize(
            vmin=self.h2.matrix.min(), vmax=self.h2.matrix.max())
        self.mappable1.set_norm(norm1)
        self.mappable2.set_norm(norm2)
        self.cbar1.set_norm(norm1)
        self.cbar2.set_norm(norm2)
        self.cbar1.update_normal(self.mappable1)
        self.cbar2.update_normal(self.mappable2)
        plt.draw()


if __name__ == "__main__":
    fn1 = '../data/cpg-islands.hg19.chr10.bed'
    fn2 = '../data/refseq.chr10.exons.bed'
    g = HilbertGUI(fn1, fn2, genome='hg19', chrom='chr10', matrix_dim=128)
    g.plot()
    plt.show()
