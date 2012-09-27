import os
import yaml
import sys
import string
from collections import defaultdict
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.widgets import Cursor, Slider, RadioButtons
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pybedtools as pbt
from scurgen.hilbert import HilbertMatrix


def data_dir():
    """
    Returns the data directory that contains example files for tests and
    documentation.
    """
    return os.path.join(os.path.dirname(__file__), 'data')


def debug_plot(h, verbose=True, nlabels=10):
    """
    Quick plot of a HilbertBase subclass that also labels the first 10 cells
    """
    imshow_kwargs = dict(
        interpolation='nearest',
        origin='upper',
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


def plot_hilbert(filenames, genome, chrom, dim=128):
    """
    Example function for a pre-configured GUI
    """
    cm_list = ['Blues', 'Reds', 'Greens', 'Greys', 'YlOrBr']
    config = dict(dim=dim, genome=genome, chrom=chrom, data=[])
    for fn, cm in zip(filenames, self.cm_list):
        config['data'].append(dict(filename=fn, colormap=cm))
    g = HilbertGUI(config)
    g.plot()
    return g


class HilbertGUI(object):
    def __init__(self, config, debug=False):
        """
        :param config:
            If a string, then treat it as a filename of a YAML config file; if
            a dictionary then treat it as the config dictionary itself.

            For each dictionary in `config['data']`, a new matrix, colorbar,
            and slider will be created using the filename and colormap
            specified.  The matrices for the files will be plotted on the same
            Axes.

            There is no limit, but colors get complicated quickly
            with, say, >3 files.

            Example config dict::

                {
                 'dim': 128,
                 'genome': 'hg19',
                 'chrom': 'chr10',
                 'data': [
                       {'filename': '../data/cpg-islands.hg19.chr10.bed',
                        'colormap': 'Blues'},

                       {'filename': '../data/refseq.chr10.exons.bed',
                        'colormap': 'Reds'}

                         ]
                }

            Example YAML file::

                dim: 128
                chrom: chr10
                genome: hg19
                data:
                    -
                        filename: ../data/cpg-islands.hg19.chr10.bed
                        colormap: Blues

                    -
                        filename: ../data/refseq.chr10.exons.bed
                        colormap: Reds


        :param debug:
            If True, then print some extra debugging info

        :param kwargs:
            Additional keyword arguments are passed to HilbertMatrix (e.g.,
            m_dim, genome, chrom)
        """
        self.config = self._parse_config(config)
        self.matrix_dim = self.config['dim']

        hilbert_matrix_kwargs = dict(
            matrix_dim=self.config['dim'],
            genome=self.config['genome'],
            chrom=self.config['chrom'])

        # self.hilberts is keyed first by chrom, then by filename; the final
        # leaves are HilbertMatrix objects
        #
        # self.hilberts = {
        #   chrom1: {
        #               filename1: HM,
        #               filename2: HM,
        #               filename3: HM,
        #           },
        #   chrom2: {
        #               filename1: HM,
        #               filename2: HM,
        #               filename3: HM,
        #           },
        # }
        #
        #
        self.hilberts = defaultdict(dict)

        # colormaps are consistent across all chroms, so it's just keyed by
        # filename:
        #
        # self.colormaps = {
        #   filename1: cmap1,
        #   filename2: cmap2,
        #   filename3: cmap3
        # }
        self.colormaps = {}

        chroms = self.config['chrom']

        if chroms == 'genome':
            chroms = pbt.chromsizes(self.config['genome']).default.keys()

        if isinstance(chroms, basestring):
            chroms = [chroms]

        for chunk in self.config['data']:
            fn = chunk['filename']
            self.colormaps[fn] = getattr(matplotlib.cm, chunk['colormap'])
            for chrom in chroms:
                hm = HilbertMatrix(fn, chrom=chrom, **hilbert_matrix_kwargs))
                hm.mask_low_values()
                self.hilberts[chrom][fn] = hm

        self.debug = debug
        self.n = len(self.config['data'])
        self.fig = plt.figure(figsize=(8, 8))

    def _parse_config(self, config):
        if isinstance(config, basestring):
            config = yaml.load(open(config))
        self._validate_config(config)
        return config

    def _validate_config(self, config):
        # TODO: more work on validation
        assert 'data' in config

    def _make_main_axes(self):
        self.ax = plt.Axes(self.fig, (0.1, 0.1, 0.8, 0.8))

        # formatting: remove ticklabels on the main axes since they're not
        # meaningful
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.divider = make_axes_locatable(self.ax)
        self.fig.add_axes(self.ax)

    def _make_colorbar_axes(self):
        cax_total_width = 0.4
        cax_total_padding = 0.8
        width = cax_total_width / self.n
        pad = cax_total_padding / self.n
        self.caxes = []
        for i in range(self.n):
            self.caxes.append(
                self.divider.append_axes('right', size=width, pad=pad))

    def _make_alpha_slider_axes(self):
        # Alpha sliders.
        self.slider_axes = []
        for i in range(self.n):
            self.slider_axes.append(
                self.divider.append_axes('bottom', size=0.1, pad=0.1))

    def _make_min_slider_axes(self):
        return
        self.min_slider_ax1 = plt.Axes(self.fig, (0.7, 0.07, 0.07, 0.02))
        self.min_slider_ax2 = plt.Axes(self.fig, (0.7, 0.02, 0.07, 0.02))
        self.fig.add_axes(self.min_slider_ax1)
        self.fig.add_axes(self.min_slider_ax2)

    def _make_annotation_axes(self):
        self.annotation_ax = plt.Axes(
            self.fig, (0.1, 0.9, 0.8, 0.05), frame_on=False)
        self.annotation_ax.set_xticks([])
        self.annotation_ax.set_yticks([])

        # necessary to capture the canvas before drawing new things to it.
        # otherwise, old chrom coordinates are left on the canvas.
        #http://stackoverflow.com/questions/6286731/ \
        #animating-matplotlib-panel-blit-leaves-old-frames-behind
        self.fig.canvas.draw()
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
        self.mappables = []

        for i in range(self.n):
            h = self.hilberts[i]
            cmap = self.colormaps[i]
            if i == 0:
                picker = 5
            else:
                picker = None
            self.mappables.append(
                self.ax.imshow(
                    h.masked, interpolation='nearest', origin='upper',
                    cmap=cmap, picker=picker))

        # Initialize alphas
        for m in self.mappables:
            m.set_alpha(0.5)

    def _matrix_colorbars(self):
        # colorbars
        self.cbars = []
        for i in range(self.n):
            m = self.mappables[i]
            cax = self.caxes[i]
            self.cbars.append(plt.colorbar(m, cax=cax))

        # Tweak colorbar labels
        for cbar in self.cbars:
            for txt in cbar.ax.get_yticklabels():
                txt.set_size(8)

    def _init_alpha_sliders(self):
        # Set up sliders with sensible default labels
        self.sliders = []

        for i in range(self.n):
            fn = self.config['data'][i]['filename']
            label = '%s: %s' % (string.letters[i], os.path.basename(fn))
            slider = Slider(
                self.slider_axes[i],
                label,
                valmin=0,
                valmax=1,
                valinit=0.5)
            slider.label.set_size(10)
            self.sliders.append(slider)

        #self.slider1.poly.set_color('#7a0510')
        #self.slider2.poly.set_color('#08316d')

    def _init_min_sliders(self):
        return
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
        for i in range(self.n):
            self.sliders[i].on_changed(
                self._slider_callback_factory(
                    self.mappables[i], self.cbars[i]))

        """
        # Min sliders change thresh
        self.min_slider1.on_changed(
            self._min_slider_callback_factory(
                self.h1, self.mappable1, self.cbar1))
        self.min_slider2.on_changed(
            self._min_slider_callback_factory(
                self.h2, self.mappable2, self.cbar2))
        """

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
            if (xi >= self.matrix_dim) or (yi >= self.matrix_dim):
                s = ""
            else:
                # Genomic coords from (x,y)
                s = '%s:%s-%s' % self.hilberts[0].xy2chrom(xi, yi)

                v = []
                for letter, h in zip(string.letters, self.hilberts):
                    v.append('%s=%s' % (letter, h.matrix[yi, xi]))
                s += self._xy_to_value_string(xi, yi)

        # Update text, redraw just the text object, and blit the background
        # previously saved
        self.current_position_label.set_text(s)
        self.annotation_ax.draw_artist(self.current_position_label)
        self.fig.canvas.blit(self.annotation_ax.bbox)

    def _xy_to_value_string(self, x, y):
        v = []
        for letter, h in zip(string.letters, self.hilberts):
            v.append('%s=%s' % (letter, h.matrix[y, x]))
        return '; '.join(v)

    def _coord_callback(self, event):
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata
        xi = int(round(x))
        yi = int(round(y))
        if self.debug:
            print
            print 'mouse x:', x, 'xi:', xi
            print 'mouse y:', y, 'yi:', yi
        s = '%s:%s-%s' % self.hilberts[0].xy2chrom(xi, yi)

        s += '\n' + self._xy_to_value_string(xi, yi)
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
        for i in range(self.n):
            norm = matplotlib.colors.LogNorm(
                vmin=self.hilberts[i].masked.min(),
                vmax=self.hilberts[i].masked.max())
            self.mappables[i].set_norm(norm)
            self.cbars[i].set_norm(norm)
            self.cbars[i].update_normal(self.mappables[i])
        plt.draw()

    def _linear(self):
        """
        Update colormaps of the plotted images to use linear-scaled color
        """
        for i in range(self.n):
            norm = matplotlib.colors.Normalize(
                vmin=self.hilberts[i].masked.min(),
                vmax=self.hilberts[i].masked.max())
            self.mappables[i].set_norm(norm)
            self.cbars[i].set_norm(norm)
            self.cbars[i].update_normal(self.mappables[i])
        plt.draw()


def gui_main(parser, args):
    g = HilbertGUI(args.config_file)
    g.plot()
    plt.show()
