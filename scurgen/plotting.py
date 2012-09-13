import os
import sys
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.widgets import Cursor, Slider, RadioButtons
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scurgen.hilbert import HilbertMatrix


def debug_plot(h, verbose=True, nlabels=10):
    """
    Quick plot of a HilbertBase subclass that also labels the first 10 cells
    """
    imshow_kwargs = dict(
        interpolation='nearest',
        origin='lower',
        cmap=matplotlib.cm.Spectral_r)
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111)
    mappable = ax.imshow(h.matrix, **imshow_kwargs)
    plt.colorbar(mappable)
    x, y, labels = h.curve()
    row, col = y, x
    ax.plot(x, y, '0.5')
    if verbose:
        for i in range(nlabels):
            # add chrom coords if it's a HilbertMatrix
            if isinstance(h, HilbertMatrix):
                label = 'i=%s\nx=%s, y=%s\nr=%s, c=%s\n%s' \
                        % (labels[i], x[i], y[i], row[i], col[i],
                           '%s:%s-%s' % h.xy2chrom(x[i], y[i]))
            else:
                label = 'i=%s\nx=%s, y=%s\nr=%s, c=%s' \
                        % (labels[i], x[i], y[i], row[i], col[i])
            ax.text(x[i], y[i], label, size=8, verticalalignment='center',
                    horizontalalignment='center')

    ax.axis('tight')
    plt.show()


class HilbertGUI(object):
    def __init__(self, intervals1, intervals2, debug=False, **kwargs):
        """
        :param intervals1:
            A file supported by pybedtools (BED, VCF, GTF, etc)

        :param intervals2:
            Another file supporte by pybedtools

        :param debug:
            If True, then print some extra debugging info

        :param kwargs:
            Additional keyword arguments are passed to HilbertMatrix (e.g.,
            m_dim, genome, chrom)
        """
        # TODO: lots of possible configuration here.  Possibly use a YAML
        # config file strategy?
        self.matrix_dim = kwargs['matrix_dim']

        self.h1 = HilbertMatrix(intervals1, **kwargs)
        self.h2 = HilbertMatrix(intervals2, **kwargs)

        self.h1.mask_low_values()
        self.h2.mask_low_values()

        self.debug = debug

        self.fig = plt.figure(figsize=(8, 8))

    def _make_main_axes(self):
        self.ax = plt.Axes(self.fig, (0.1, 0.1, 0.8, 0.8))

        # formatting: remove ticklabels on the main axes since they're not
        # meaningful
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.fig.add_axes(self.ax)

    def _make_colorbar_axes(self):
        # from axes_grid toolkit
        # TODO: eventually support n-way comparisons
        divider = make_axes_locatable(self.ax)
        self.cax1 = divider.append_axes('right', size=0.2, pad=0.3)
        self.cax2 = divider.append_axes('right', size=0.2, pad=0.3)

    def _make_alpha_slider_axes(self):
        # Alpha sliders.
        # TODO: eventually n-way comparisons; these should be appended on the
        # bottom of self.ax?
        self.slider_ax1 = plt.Axes(self.fig, (0.3, 0.07, 0.3, 0.02))
        self.slider_ax2 = plt.Axes(self.fig, (0.3, 0.02, 0.3, 0.02))
        self.fig.add_axes(self.slider_ax1)
        self.fig.add_axes(self.slider_ax2)

    def _make_min_slider_axes(self):
        self.min_slider_ax1 = plt.Axes(self.fig, (0.7, 0.07, 0.07, 0.02))
        self.min_slider_ax2 = plt.Axes(self.fig, (0.7, 0.02, 0.07, 0.02))
        self.fig.add_axes(self.min_slider_ax1)
        self.fig.add_axes(self.min_slider_ax2)

    def _make_annotation_axes(self):
        self.annotation_ax = plt.Axes(
            self.fig, (0.1, 0.9, 0.8, 0.05), frame_on=False)
        self.annotation_ax.set_xticks([])
        self.annotation_ax.set_yticks([])
        self.background = self.fig.canvas.copy_from_bbox(
            self.annotation_ax.bbox)
        self.current_position_label = self.annotation_ax.text(
            .5, .5, 'position...', horizontalalignment='center',
            verticalalignment='center', size=10, animated=True)
        self.fig.add_axes(self.annotation_ax)

    def _make_radio_axes(self):
        self.radio_ax = plt.Axes(self.fig, (0.85, 0.02, 0.1, 0.1))
        self.fig.add_axes(self.radio_ax)

    def _imshow_matrices(self):
        # plot the matrices on top of each other. Note that only one gets the
        # `picker` kwarg; otherwise the callback will trigger multiple times.
        self.mappable1 = self.ax.imshow(
            self.h1.masked, interpolation='nearest',
            origin='lower',
            cmap=matplotlib.cm.Reds, picker=5)

        self.mappable2 = self.ax.imshow(
            self.h2.masked, interpolation='nearest',
            origin='lower',
            cmap=matplotlib.cm.Blues)

        # Initialize alphas
        self.mappable1.set_alpha(0.5)
        self.mappable2.set_alpha(0.5)

    def _matrix_colorbars(self):
        # colorbars
        self.cbar1 = plt.colorbar(self.mappable1, cax=self.cax1)
        self.cbar2 = plt.colorbar(self.mappable2, cax=self.cax2)

        # Tweak colorbar labels
        for cbar in [self.cbar1, self.cbar2]:
            for txt in cbar.ax.get_yticklabels():
                txt.set_size(8)

    def _init_alpha_sliders(self):
        # Set up sliders with sensible default labels
        self.slider1 = Slider(
            self.slider_ax1,
            'a: ' + os.path.basename(self.h1.file),
            valmin=0,
            valmax=1,
            valinit=0.5)

        self.slider2 = Slider(
            self.slider_ax2,
            'b: ' + os.path.basename(self.h2.file),
            valmin=0,
            valmax=1,
            valinit=0.5)

        self.slider1.poly.set_color('#7a0510')
        self.slider2.poly.set_color('#08316d')
        self.slider1.label.set_size(10)
        self.slider2.label.set_size(10)

    def _init_min_sliders(self):
        self.min_slider1 = Slider(
            self.min_slider_ax1,
            'min',
            valmin=self.h1.masked.min(),
            valmax=self.h1.masked.max(),
            valinit=0)

        self.min_slider2 = Slider(
            self.min_slider_ax2,
            'min',
            valmin=self.h2.masked.min(),
            valmax=self.h2.masked.max(),
            valinit=0)

        self.min_slider1.poly.set_color('#7a0510')
        self.min_slider2.poly.set_color('#08316d')
        self.min_slider1.label.set_size(10)
        self.min_slider2.label.set_size(10)

    def _init_radio(self):
        # For controlling log/linear color scale
        self.radio = RadioButtons(
            self.radio_ax,
            labels=['log', 'linear'],
            active=1)

    def _make_connections(self):
        # Alpha sliders
        self.slider1.on_changed(
            self._slider_callback_factory(self.mappable1, self.cbar1))
        self.slider2.on_changed(
            self._slider_callback_factory(self.mappable2, self.cbar2))

        # Min sliders change thresh
        self.min_slider1.on_changed(
            self._min_slider_callback_factory(
                self.h1, self.mappable1, self.cbar1))
        self.min_slider2.on_changed(
            self._min_slider_callback_factory(
                self.h2, self.mappable2, self.cbar2))

        # Radio callback changes color scale
        self.radio.on_clicked(self._radio_callback)

        self.fig.canvas.mpl_connect('motion_notify_event', self._coord_tracker)
        self.fig.canvas.mpl_connect('pick_event', self._coord_callback)

    def plot(self):
        """
        Does most of the work to set up the figure, axes, and widgets.
        """
        # These methods construct axes in the right places
        self._make_main_axes()
        self._make_colorbar_axes()
        self._make_alpha_slider_axes()
        self._make_min_slider_axes()
        self._make_annotation_axes()
        self._make_radio_axes()

        # Plot the matrices and their colorbars
        self._imshow_matrices()
        self._matrix_colorbars()

        # Initialize the various widgets
        self._init_alpha_sliders()
        self._init_min_sliders()
        self._init_radio()

        # Connect callbacks to events
        self._make_connections()

    def _coord_tracker(self, event):
        """
        Callback that updates text based on the genomic coords of the current
        mouse position.
        """
        # Restore original background
        self.fig.canvas.restore_region(self.background)

        # These will be None if the mouse is not in an Axes, so the string
        # should be empty.
        #
        # Also, make sure we're in the imshow Axes -- don't want crazy genomic
        # coords from being in the colorbar or annotation Axes
        x = event.xdata
        y = event.ydata
        if (x is None) or (y is None) or (event.inaxes is not self.ax):
            s = ""

        # Get matrix coords
        else:
            xi = int(round(x))
            yi = int(round(y))

            # If you move the mouse off the edge of the main Axes, rounding to
            # the nearest row or col may get you a value that's greater than
            # the number of rows/cols.  In this case, treat it similar to being
            # out of the Axes, with an empty string.
            if (xi >= self.h1.m_dim) or (yi >= self.h1.m_dim):
                s = ""
            else:
                # Genomic coords from (x,y)
                s = '%s:%s-%s' % self.h1.xy2chrom(xi, yi)

                # Values of the underlying matrices.  Note xi,yi swap for row,
                # col coords
                s += ' [a=%s; b=%s]' \
                    % (self.h1.matrix[yi, xi], self.h2.matrix[yi, xi])

        # Update text, redraw just the text object, and blit the background
        # previously saved
        self.current_position_label.set_text(s)
        self.annotation_ax.draw_artist(self.current_position_label)
        self.fig.canvas.blit(self.annotation_ax.bbox)

    def _coord_callback(self, event):
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata
        xi = int(round(x))
        yi = int(round(y))
        if self.debug:
            print
            print 'mouse x:', x, 'xi:', xi
            print 'mouse y:', y, 'yi:', yi
        s = '%s:%s-%s' % self.h1.xy2chrom(xi, yi)

        # Note that xi=cols and yi=rows, hence the indexing switcharoo
        s += ' [a=%s; b=%s]' % (self.h1.matrix[yi, xi], self.h2.matrix[yi, xi])
        print s
        sys.stdout.flush()

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

    def _min_slider_callback_factory(self, h, mappable, cbar):
        def _min_slider_callback(x):
            h.mask_low_values(min_val=x)
            mappable.set_data(h.masked)

        return _min_slider_callback

    def _log(self):
        """
        Update colomaps of the plotted images to use log-scaled color
        """
        norm1 = matplotlib.colors.LogNorm(
            vmin=self.h1.masked.min(), vmax=self.h1.masked.max())
        norm2 = matplotlib.colors.LogNorm(
            vmin=self.h2.masked.min(), vmax=self.h2.masked.max())
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
            vmin=self.h1.masked.min(), vmax=self.h1.masked.max())
        norm2 = matplotlib.colors.Normalize(
            vmin=self.h2.masked.min(), vmax=self.h2.masked.max())
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
