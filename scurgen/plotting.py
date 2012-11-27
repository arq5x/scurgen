import os
import yaml
import sys
import string
from collections import defaultdict
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gs
import pybedtools as pbt
from scurgen.hilbert import HilbertMatrix, HilbertMatrixBigWig
from scipy import stats


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


class HilbertPlot(object):
    def __init__(self, config, debug=False):
        """
        Class for handling plotting multiple Hilbert matrices.

        This class is designed for programmatic access; if you want to interact
        with it via a GUI, then use the HilbertGUI subclass which adds the GUI
        elements.

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
            genome=self.config['genome'])

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

        if chroms == 'chroms':
            chroms = pbt.chromsizes(self.config['genome']).default.keys()

        if isinstance(chroms, basestring):
            chroms = [chroms]

        self.chroms = chroms
        self.fns = []
        for chunk in self.config['data']:
            is_bigwig = chunk.get('bigwig', False)
            if is_bigwig:
                HilbertClass = HilbertMatrixBigWig
            else:
                HilbertClass = HilbertMatrix

            fn = chunk['filename']
            self.fns.append(fn)
            self.colormaps[fn] = getattr(matplotlib.cm, chunk['colormap'])
            for chrom in self.chroms:
                hm = HilbertClass(fn, chrom=chrom, **hilbert_matrix_kwargs)
                hm.mask_low_values()
                self.hilberts[chrom][fn] = hm

        self.debug = debug
        self.nfiles = len(self.config['data'])
        self.nchroms = len(chroms)
        self.annotation_ax = None

    def _parse_config(self, config):
        """
        Helper method to parse the configuration.

        Returns a validated config dictionary.  `config` can be a filename of
        a YAML-formatted file or an existing dictionary
        """
        if isinstance(config, basestring):
            config = yaml.load(open(config))
        self._validate_config(config)
        return config

    def _validate_config(self, config):
        """
        Helper method to validate the configuration
        """

        # TODO: more work on validation
        assert 'data' in config

    def _configure_axes(self):
        """
        Helper method to construct axes.

        Given the number of chromosomes and the number of data files
        configured, this creates a bunch of subplots.

        Abuses matplotlib.gridspec.GridSpec
        """
        # The area that will be subdivided into subplots per chrom.
        # TODO: expose this to the figure-saving method, which should get its
        # bbox [in part] from these coords.
        CHROM = dict(
            left=0.05,
            right=0.8,
            top=0.9,
            bottom=0.2,
            wspace=0.01,
            hspace=0.1)

        # Area used for colorbars
        # TODO: similar to CHROM, expose this to figure-saving method
        CBAR_PAD = 0.01
        CBAR = dict(
            left=CHROM['right'] + CBAR_PAD,
            right=0.95,
            wspace=1.5,
            top=CHROM['top'],
            bottom=CHROM['bottom'],
        )

        # Set up the grids upon which new axes will eventually be organized and
        # created
        nrows = int(np.ceil(np.sqrt(self.nchroms)))
        ncols = nrows

        assert nrows * ncols >= len(self.chroms)

        chrom_grid = gs.GridSpec(nrows, ncols)
        chrom_grid.update(**CHROM)

        colorbar_grid = gs.GridSpec(1, self.nfiles)
        colorbar_grid.update(**CBAR)

        # Now create the actual axes, and store them in self.*_axs dictionaries
        # for later access -- keyed by chrom (for chrom_axs) or filenames (all
        # others)
        self.chrom_axs = {}
        for chrom, spec in zip(self.chroms, chrom_grid):
            ax = plt.subplot(spec)
            ax.set_yticks([])
            ax.set_xticks([])
            ax.set_title(chrom, size=10)
            self.chrom_axs[chrom] = ax

        self.colorbar_axs = {}
        for fn, spec in zip(self.fns, colorbar_grid):
            self.colorbar_axs[fn] = plt.subplot(spec)

    def _make_annotation_axes(self, event=None):
        """
        Helper method to make an axes upon which to print genomic coords
        """
        self.annotation_ax = plt.Axes(
            self.fig, (0.1, 0.95, 0.8, 0.03), frame_on=False)
        self.annotation_ax.set_xticks([])
        self.annotation_ax.set_yticks([])

        # necessary to capture the canvas before drawing new things to it.
        # otherwise, old chrom coordinates are left on the canvas.
        #http://stackoverflow.com/questions/6286731/ \
        #animating-matplotlib-panel-blit-leaves-old-frames-behind
        self.fig.canvas.draw()
        self._background_snapshot()
        self.current_position_label = self.annotation_ax.text(
            .5, .5, 'position...', horizontalalignment='center',
            verticalalignment='center', size=10, animated=True,
            transform=self.annotation_ax.transAxes)
        self.fig.add_axes(self.annotation_ax)

    def _background_snapshot(self, event=None):
        """
        Helper method to take a snapshot of the annotation axes
        """
        self.background = self.fig.canvas.copy_from_bbox(
            self.annotation_ax.bbox)

    def _imshow_matrices(self):
        """
        Helper method to plot matrices.

        Assumes HilbertMatrix objects have already been created in
        self.hilberts; this just imshow()s each underlying matrix on the
        appropriate axes
        """
        # Like the structure of self.hilberts, self.mappables is keyed by chrom
        # and then filename.
        self.mappables = defaultdict(dict)

        for chrom in self.chroms:
            for fn in self.fns:
                h = self.hilberts[chrom][fn]
                if len(self.mappables[chrom]) == 0:
                    picker = 5
                else:
                    picker = None
                cmap = self.colormaps[fn]
                ax = self.chrom_axs[chrom]
                mappable = ax.imshow(
                    h.masked, interpolation='nearest', origin='upper',
                    cmap=cmap, picker=picker)
                mappable.set_alpha(0.5)
                self.mappables[chrom][fn] = mappable

    def _matrix_colorbars(self):
        """
        Helper method to add colorbars.

        Assumes that self._configure_axes() and self._imshow_matrices() have
        been called so that self.colorbar_axs and self.mappables have been
        populated.
        """
        self.colorbars = {}

        # even though we have nchroms x nfiles mappables, we only need to make
        # colorbars for nfiles of them.  So just grab the mappables for the
        # first configured chrom.
        chrom = self.chroms[0]
        for fn in self.fns:
            self.colorbars[fn] = plt.colorbar(
                self.mappables[chrom][fn],
                cax=self.colorbar_axs[fn])

        # Tweak colorbar labels
        for cbar in self.colorbars.values():
            for txt in cbar.ax.get_yticklabels():
                txt.set_size(8)

    def plot(self, figsize=(8, 6)):
        """
        Plots and shows the configured HilbertMatrix objects.

        This method calls all of the plotting setup helper methods which in
        turn set up all the axes needed for plotting.

        :param figsize:
            (width, height) tuple for the figure size to create
        """
        # These methods construct axes in the right places
        self.fig = plt.figure(figsize=figsize)
        self._configure_axes()
        self._make_annotation_axes()

        # Plot the matrices and their colorbars
        self._imshow_matrices()
        self._matrix_colorbars()
        self.fig.canvas.mpl_connect('motion_notify_event', self._coord_tracker)
        self.fig.canvas.mpl_connect('pick_event', self._coord_callback)
        self.fig.canvas.mpl_connect('resize_event', self._background_snapshot)

    # Helper methods for getting various info when given an axes --------------

    def _chrom_from_axes(self, ax):
        """
        Helper method to return which chrom an axes represents

        Reverse lookup into self.chrom_axs to return what chromsome a given
        axes represents.
        """
        for chrom, chrom_ax in self.chrom_axs.items():
            if chrom_ax == ax:
                return chrom

    def _hilberts_from_axes(self, ax):
        """
        Helper method to return HilbertMatrix objects plotted on an axes.

        Given an axes, return the {filename: HilbertMatrix} dictionary of
        HilbertMatrix objects that are plotted on it.
        """
        chrom = self._chrom_from_axes(ax)
        return self.hilberts[chrom]

    def _first_hilbert_from_axes(self, ax):
        """
        Helper method to return a single HilbertMatrix

        Some of the callbacks don't care which of the possibly multiple
        HilbertMatrix objects are plotted on the matrix -- they just need one
        of them in order to get the genomic coords.
        """
        return self._hilberts_from_axes(ax).itervalues().next()

    # Callbacks and callback helpers ------------------------------------------

    def _coord_tracker(self, event):
        """
        Callback to show genomic coords of the current mouse position.
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
        if (x is None) or (y is None) \
                or (event.inaxes not in self.chrom_axs.values()):
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
                # Identify the HilbertMatrix objects for this axes
                hms = self._hilberts_from_axes(event.inaxes)

                # Really we only need one of them.
                hm = hms.itervalues().next()

                # Genomic coords from (x,y)
                s = '%s:%s-%s' % hm.xy2chrom(xi, yi)

        # Update text, redraw just the text object, and blit the background
        # previously saved
        self.current_position_label.set_text(s)
        self.annotation_ax.draw_artist(self.current_position_label)
        self.fig.canvas.blit(self.annotation_ax.bbox)

    def _coord_callback(self, event):
        """
        Callback function to report chrom/start/stop data of the mouse pointer
        """
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata
        xi = int(round(x))
        yi = int(round(y))
        if self.debug:
            print
            print 'mouse x:', x, 'xi:', xi
            print 'mouse y:', y, 'yi:', yi
        h = self._first_hilbert_from_axes(event.artist.axes)
        s = '%s:%s-%s' % h.xy2chrom(xi, yi)
        print s
        sys.stdout.flush()

    def set_alpha(self, fn, alpha):
        """
        Set alpha for file `fn`

        :param fn:
            Which filename to set the alpha for

        :param alpha:
            Transparency of the matrix colors (0-1)
        """
        if isinstance(fn, int):
            fn = self.fns[fn]
        for chrom in self.chroms:
            mappable = self.mappables[chrom][fn]
            cbar = self.colorbars[fn]
            mappable.set_alpha(alpha)
            cbar.set_alpha(alpha)
            cbar.update_normal(mappable)
        plt.draw()

    def _colormap_normalizer(self, fn, norm):
        """
        Helper function to apply the norm to all matrices for filename.
        """
        for chrom in self.chroms:
            self.mappables[chrom][fn].set_norm(norm)
        self.colorbars[fn].set_norm(norm)
        #self.colorbars[fn].update_normal(self.mappables[chrom][fn])

    def _min_max_for_fn(self, fn):
        """
        Helper function to get data limits.

        Identifies the (masked) min/max values across all chroms for filename
        `fn`
        """
        mns, mxs = [], []
        for chrom in self.chroms:
            m = self.hilberts[chrom][fn].masked
            mns.append(m.min())
            mxs.append(m.max())
        mn = min(mns)
        mx = max(mxs)
        return mn, mx

    def set_log(self):
        """
        Update colomaps of the plotted images to use log-scaled color
        """
        for fn in self.fns:
            mn, mx = self._min_max_for_fn(fn)
            norm = matplotlib.colors.LogNorm(vmin=mn, vmax=mx)
            self._colormap_normalizer(fn, norm)
        plt.draw()

    def set_linear(self):
        """
        Update colormaps of the plotted images to use linear-scaled color
        """
        for fn in self.fns:
            mn, mx = self._min_max_for_fn(fn)
            norm = matplotlib.colors.Normalize(vmin=mn, vmax=mx)
            self._colormap_normalizer(fn, norm)
        plt.draw()

    def set_colorbar_limits(self, min_val, max_val, fn=None, percentile=False):
        """
        Manually set colorbar limits.

        :param min_val: Minimum value
        :param max_val: Maximum value
        :param fn:
            If not None, then apply the colorbar limits only to the matrices
            showing this filename; otherwise apply to all files shown.
        :param percentile:
            If True, then assume `min_val` and `max_val` represent percentiles
            (e.g., 95 for 95th percentile) and re-calculate actual values based
            on the data
        """
        if fn is not None:
            fns = [fn]
        else:
            fns = self.fns
        for fn in fns:
            to_concat = []
            for chrom in self.chroms:
                if percentile:
                    to_concat.append(self.hilberts[chrom][fn].matrix)

            if percentile:
                concatenated = np.concatenate([i.ravel() for i in to_concat])
                vmin = stats.scoreatpercentile(concatenated, min_val)
                vmax = stats.scoreatpercentile(concatenated, max_val)

            else:
                vmin = min_val
                vmax = max_val
            norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
            self._colormap_normalizer(fn, norm)

        plt.draw()

    def reset_mask(self, fn=None):
        """
        Reset the mask

        :param fn:
            If not None, then reset the mask only for the matrices showing this
            filename; otherwise apply to all matrices.
        """
        if fn is not None:
            fns = [fn]
        else:
            fns = self.fns

        for chrom in self.chroms:
            for fn in fns:
                h = self.hilberts[chrom][fn]
                h.masked.mask[:] = False
                mn, mx = self._min_max_for_fn(fn)
                norm = matplotlib.colors.Normalize(vmin=mn, vmax=mx)
                self._colormap_normalizer(fn, norm)
        plt.draw()

    def mask_low_values(self, min_val=0, fn=None):
        """
        Mask low values across all matrices

        Colorbar limits will be re-calculated to reflect the newly masked data.

        :param min_val:
            Values <= `min_val` will be masked out.

        :param fn:
            If not None, then only apply the masking to filename `fn` matrices
            (across all chromsomes)
        """
        if fn is not None:
            fns = [fn]
        else:
            fns = self.fns
        for chrom in self.chroms:
            for fn in fns:
                h = self.hilberts[chrom][fn]
                h.mask_low_values(min_val)
                mn, mx = self._min_max_for_fn(fn)
                norm = matplotlib.colors.Normalize(vmin=mn, vmax=mx)
                self._colormap_normalizer(fn, norm)
        plt.draw()

    def mask_high_values(self, max_val=np.inf):
        """
        Mask high values across all matrices

        Colorbar limits will be re-calculated to reflect the newly masked data.

        :param max_val:
            Values >= `max_val` will be masked out.

        :param fn:
            If not None, then only apply the masking to filename `fn` matrices
            (across all chromsomes)
        """
        if fn is not None:
            fns = [fn]
        else:
            fns = self.fns
        for chrom in self.chroms:
            for fn in fns:
                h = self.hilberts[chrom][fn]
                h.mask_high_values(max_val)
                mn, mx = self._min_max_for_fn(fn)
                norm = matplotlib.colors.Normalize(vmin=mn, vmax=mx)
                self._colormap_normalizer(fn, norm)
        plt.draw()

    def save(self, *args, **kwargs):
        """
        Save the figure.

        Computes the bounding box of the chrom axes and colorbar axes alone, so
        that any additional GUI elements are not saved.

        :param fname:
            File name to save as; extension is meaningful

        Additional args and kwargs are passed verbatim to
        matplotlib.pyplot.savefig
        """
        # Figure.savefig() has a `bbox_inches` kwarg.  The trick is to figure
        # out, given the *display* coordinates covered by the chrom axes and
        # colorbar axes, what the current figure *inches* should be.

        # This gets the union of all bboxes of axes we care about
        artists = []
        axes = self.chrom_axs.values() + self.colorbar_axs.values()
        bboxes = []
        for ax in axes:
            bboxes.append(ax.get_window_extent())
            for txt in ax.findobj(matplotlib.text.Text):
                bboxes.append(
                    txt.get_window_extent(
                        renderer=self.fig.canvas.renderer)
                )

        union = matplotlib.transforms.Bbox.union(bboxes)

        # Convert to inches, thanks to
        # http://stackoverflow.com/questions/4325733/\
        # save-a-subplot-in-matplotlib
        extent = union.transformed(self.fig.dpi_scale_trans.inverted())
        self.fig.savefig(bbox_inches=extent, *args, **kwargs)


class HilbertGUI(HilbertPlot):
    def __init__(self, *args, **kwargs):
        """
        Similar to HilbertPlot class, but automatically adds GUI elements.
        """
        super(HilbertGUI, self).__init__(*args, **kwargs)

    def _configure_axes(self):
        """
        Given the number of chromosomes and the number of data files
        configured, create a bunch of subplots.

        Abuses matplotlib.gridspec.GridSpec
        """
        # The area that will be subdivided into subplots per chrom.
        # TODO: expose this to the figure-saving method, which should get its
        # bbox [in part] from these coords.
        CHROM = dict(
            left=0.05,
            right=0.8,
            top=0.9,
            bottom=0.2,
            wspace=0.01,
            hspace=0.1)

        # Area used for alpha sliders
        SLIDER_PAD = 0.01
        SLIDER = dict(
            left=0.15,
            right=0.70,
            bottom=0.1,
            top=CHROM['bottom'] - SLIDER_PAD,
            hspace=0.5,
        )

        # Area used for colorbars
        # TODO: similar to CHROM, expose this to figure-saving method
        CBAR_PAD = 0.01
        CBAR = dict(
            left=CHROM['right'] + CBAR_PAD,
            right=0.95,
            wspace=1.5,
            top=CHROM['top'],
            bottom=CHROM['bottom'],
        )

        # Area used for checkboxes.
        CHECKS = dict(
            top=SLIDER['top'],
            bottom=SLIDER['bottom'],
            left=CBAR['left'],
            right=CBAR['right'],
            wspace=CBAR['wspace'],
            hspace=SLIDER['hspace'])

        RADIO_PAD = 0.01
        RADIO = (
            CHROM['left'],    # left
            SLIDER['bottom'],  # bottom
            SLIDER['left'] - CHROM['left'] - RADIO_PAD,  # width
            SLIDER['top'] - SLIDER['bottom'],    # height
        )

        self.radio_ax = plt.Axes(self.fig, RADIO)
        self.fig.add_axes(self.radio_ax)

        # Set up the grids upon which new axes will eventually be organized and
        # created
        nrows = int(np.ceil(np.sqrt(self.nchroms)))
        ncols = nrows

        assert nrows * ncols >= len(self.chroms)

        chrom_grid = gs.GridSpec(nrows, ncols)
        chrom_grid.update(**CHROM)

        slider_grid = gs.GridSpec(self.nfiles, 1)
        slider_grid.update(**SLIDER)

        colorbar_grid = gs.GridSpec(1, self.nfiles)
        colorbar_grid.update(**CBAR)

        checks_grid = gs.GridSpec(self.nfiles, self.nfiles)
        checks_grid.update(**CHECKS)

        # Now create the actual axes, and store them in self.*_axs dictionaries
        # for later access -- keyed by chrom (for chrom_axs) or filenames (all
        # others)
        self.chrom_axs = {}
        for chrom, spec in zip(self.chroms, chrom_grid):
            ax = plt.subplot(spec)
            ax.set_yticks([])
            ax.set_xticks([])
            ax.set_title(chrom, size=10)
            self.chrom_axs[chrom] = ax

        self.slider_axs = {}
        for fn, spec in zip(self.fns, slider_grid):
            self.slider_axs[fn] = plt.subplot(spec)

        self.colorbar_axs = {}
        for fn, spec in zip(self.fns, colorbar_grid):
            self.colorbar_axs[fn] = plt.subplot(spec)

        self.checks_axs = {}
        for i in range(self.nfiles):
            # this way only axes on the diagonal across the checkbox grid will
            # be created
            ax = plt.subplot(checks_grid[i, i])
            ax.set_yticks([])
            ax.set_xticks([])
            ax.patch.set_color('k')
            ax.SHOW = True
            ax.last_slider_val = 0.5
            self.checks_axs[self.fns[i]] = ax

    # Initialize widgets ------------------------------------------------------

    def _init_alpha_sliders(self):
        """
        Add sliders to self.slider_axs.  Assumes self._configure_axes has been
        called so that self.slider_axs has been populated.
        """
        self.sliders = {}

        for fn in self.fns:
            # Disable for now....
            #label = '%s: %s' % (string.letters[i], os.path.basename(fn))
            label = ""
            slider = Slider(
                self.slider_axs[fn],
                label,
                valmin=0,
                valmax=1,
                valinit=0.5,
                facecolor='0.3')
            slider.label.set_size(10)
            self.sliders[fn] = slider

        # TODO: it would be a nice touch if slider color == max cmap color

    def _init_radio(self):
        """
        Add radio buttons to the radio ax for color scale.
        self._make_radio_axes() needs to have been called.
        """
        self.radio = RadioButtons(
            self.radio_ax,
            labels=['log', 'linear'],
            active=1)

    def _make_connections(self):
        # Alpha sliders
        for fn in self.fns:
            self.sliders[fn].on_changed(
                self._slider_callback_factory(fn))

        # Radio callback changes color scale
        self.radio.on_clicked(self._radio_callback)

        self.fig.canvas.mpl_connect('button_press_event', self._check_callback)
        self.fig.canvas.mpl_connect('motion_notify_event', self._coord_tracker)
        self.fig.canvas.mpl_connect('pick_event', self._coord_callback)

    # Plot data ---------------------------------------------------------------

    def plot(self, figsize=(8, 8)):
        """
        Does most of the work to set up the figure, axes, and widgets.
        """
        super(HilbertGUI, self).plot(figsize=figsize)

        # Initialize the various widgets
        self._init_alpha_sliders()
        self._init_radio()

        # Connect callbacks to events
        self._make_connections()

    # Helper methods for getting various info when given an axes --------------

    def _radio_callback(self, label):
        # Hello? Am I the 9th caller!?
        if label == 'log':
            self.set_log()
        elif label == 'linear':
            self.set_linear()
        else:
            raise ValueError("unspecified label for radio button")

    def _check_callback(self, event):
        """
        Toggles the alpha between 0.0 and 0.5 and also changes the color of the
        checkbox axes.
        """
        ax = event.inaxes
        fn = self._fn_from_check_ax(ax)
        if ax in self.checks_axs.values():
            ax.SHOW = not ax.SHOW
            color = {True: 'k', False: 'w'}
            if ax.SHOW:
                self.sliders[fn].set_val(ax.last_slider_val)
            else:
                ax.last_slider_val = self.sliders[fn].val
                self.sliders[fn].set_val(0)
            ax.patch.set_color(color[ax.SHOW])

            # explicitly update all the matrices with the new value
            for chrom in self.chroms:
                self.mappables[chrom][fn].set_alpha(self.sliders[fn].val)

        plt.draw()

    def _fn_from_check_ax(self, ax):
        for fn, check_ax in self.checks_axs.iteritems():
            if check_ax == ax:
                return fn

    def _slider_callback_factory(self, fn):
        """
        Given a filename, return a function that will modify all that
        filename's mappables based on the slider's value.
        """
        def _slider_callback(x):
            self.set_alpha(fn, x)

        return _slider_callback

    def set_alpha(self, fn, alpha):
        """
        Set alpha for file `fn`

        :param fn:
            Which filename to set the alpha for

        :param alpha:
            Transparency of the matrix colors (0-1)
        """
        super(HilbertGUI, self).set_alpha(fn, alpha)
        self.sliders[fn].val = alpha


def gui_main(parser, args):
    g = HilbertGUI(args.config_file)
    g.plot()
    plt.show()


def _debug_HilbertPlot():
    config = dict(
        dim=16,
        chrom=['chr10', 'chr11'],
        genome='hg19',
        data=[
            dict(
                filename='data/cpg-islands.hg19.chr10.bed',
                colormap='Blues'),
            dict(
                filename='data/refseq.chr10.exons.bed',
                colormap='Reds'),
            dict(
                filename='data/phastcons.chr10.bed',
                colormap='Spectral_r')
        ]
    )
    g = HilbertPlot(config)
    g.plot()
    plt.show()
    return g


def _debug_HilbertGUI():
    config = dict(
        dim=16,
        chrom=['chr10', 'chr11'],
        genome='hg19',
        data=[
            dict(
                filename='data/cpg-islands.hg19.chr10.bed',
                colormap='Blues'),
            dict(
                filename='data/refseq.chr10.exons.bed',
                colormap='Reds'),
            dict(
                filename='data/phastcons.chr10.bed',
                colormap='Spectral_r')
        ]
    )

    g = HilbertGUI(config)
    g.plot()
    plt.show()
    return g

if __name__ == '__main__':
    import sys
    try:
        to_run = sys.argv[1]
        d = {'Plot': _debug_HilbertPlot,
             'GUI': _debug_HilbertGUI}

        g = d[to_run]()
    except IndexError:
        print 'first arg is Plot or GUI'
