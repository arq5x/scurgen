from bx.bbi.bigwig_file import BigWigFile
import numpy as np
import pybedtools as pbt
import sys
import os
import re
import matplotlib
import subprocess
import bisect
from pybedtools import genome_registry


def rot(n, x, y, rx, ry):
    if (ry == 0):
        if (rx == 1):
            x = n - 1 - x
            y = n - 1 - y
        # Swap x and y
        x, y = y, x
    return (x, y)


def d2xy(n, d):
    """
    Convert a distance into (x,y) coords of the matrix of dimension `n`
    """
    t = d
    x = 0
    y = 0
    s = 1
    while s < n:
        rx = 1 & (t / 2)
        ry = 1 & (t ^ rx)
        (x, y) = rot(s, x, y, rx, ry)
        x += (s * rx)
        y += (s * ry)
        t /= 4
        s *= 2
    return (x, y)


def xy2d(n, x, y):
    """
    Convert an (x, y) coordinate into distance
    """
    d = 0
    s = n / 2
    while s > 0:
        rx = (x & s) > 0
        ry = (y & s) > 0
        d += s * s * ((3 * rx) ^ ry)
        (x, y) = rot(s, x, y, rx, ry)
        s /= 2
    return d


def rc2d(n, row, col):
    """
    Same as xy2d, but use row, column semantics instead of x, y
    """
    x = col
    y = n - row - 1
    return xy2d(n, x, y)


def d2rc(n, d):
    """
    Same as d2xy, but use row, column semantics instead of x, y
    """
    x, y = d2xy(n, d)
    col = x
    row = y
    return row, col


def get_interval_from_string(s):
    """
    e.g., ``"chr21:5-100"`` -> ``('chr2l', 5, 100)``
    """
    chrom_range_pattern = re.compile('(\S+)\:([0-9]+)\-([0-9]+)')

    if chrom_range_pattern.search(s):
        (chrom, start, end) = re.findall(chrom_range_pattern, s)[0]
        start = int(start)
        end = int(end)
        return (chrom, start, end)
    else:
        return None


class HilbertBase(object):
    def __init__(self, matrix_dim):
        """
        Initializes a basic Hilbert curve.

        :param matrix_dim:
            The number of dimensions on a side (power of 2)
        """
        if np.log2(matrix_dim) % 1 != 0:
            raise ValueError('matrix_dim %s not a power of 2' % matrix_dim)
        self.matrix_dim = matrix_dim
        self.ncells = matrix_dim * matrix_dim
        self.matrix = np.zeros(
            (self.matrix_dim, self.matrix_dim),
            dtype=np.float
        )
        self.masked = np.ma.masked_equal(self.matrix, 1)

    def update(self, d1, d2, value=1, func=np.add):
        """
        Update the matrix between cells `d1` and `d2` (inclusive) by `value`.

        :param d1:
            First cell to update

        :param d2:
            Last cell to update

        :param value:
            Value to update by

        :param func:
            Optional arbitrary function that returns a float.  It should have
            the signature::

                func(existing_value, value)

            By default, func=np.add, so the cells will be simply incremented by
            `value`.
        """
        for dist in xrange(d1, d2 + 1):
            # filling in matrix, so use row, col coords
            row, col = d2rc(self.matrix_dim, dist)
            self.matrix[row, col] = func(self.matrix[row, col], value)

    def mask_low_values(self, min_val=0):
        """
        Mask values <= `min_val`
        """
        self.masked = np.ma.masked_where(self.masked.mask | (self.matrix <= min_val), self.matrix)

    def mask_high_values(self, max_val=np.inf):
        """
        Mask values >= `min_val`
        """
        self.masked = np.ma.masked_where(self.masked.mask | (self.matrix >= max_val), self.matrix)

    def reset(self):
        """
        Resets the matrix to zeros everywhere
        """
        self.matrix[:] = 0.0

    def curve(self):
        """
        Returns a 3-tuple of the x-coords, y-coords, and labels for the curve.
        The (x, y) coords correspond to (col, row) positions in a matrix (note
        the switch) such that that the following will work::

            x, y, labels = x.curve()
            plt.imshow(x.matrix)
            plt.plot(x, y)

        In this case the (x, y) coords fall within the center of each cell.
        """
        # To be plotted, so use x, y coords instead of row, col coords
        xs, ys = [], []
        for i in range(self.ncells):
            x, y = d2xy(self.matrix_dim, i)
            xs.append(x)
            ys.append(y)

        return np.array(xs), np.array(ys), range(self.ncells)


class HilbertNormalized(HilbertBase):
    def __init__(self, matrix_dim, length):
        """
        Subclass of HilbertBase to handle distance in arbitrary units
        (e.g., genomic bp) instead of in cell numbers.

        :param matrix_dim:
            The number of dimensions on a side (power of 2)

        :param length:
            The total length represented by this curve
        """
        super(HilbertNormalized, self).__init__(matrix_dim)
        self.length = length
        self.norm_factor = self.length / float(self.ncells)

    @property
    def dist_per_cell(self):
        """
        Returns the distance represented by each cell in the matrix.

        (alias for self.norm_factor)
        """
        return self.norm_factor

    def normalize(self, d):
        """
        Convert distance from bp to number of cells.
        """
        return int(d / self.norm_factor)

    def update(self, d1, d2, value=1, func=np.add, cells=False):
        """
        Update the matrix between distances `d1` and `d2` (inclusive) by
        `value`.

        :param d1:
            Beginning of the distance to update, or, if cells=True, assume `d1`
            has already been normalized to a cell distance

        :param d2:
            End of the distance to update, or, if cells=True, assume `d2` has
            already been normalized to a cell distance

        :param value:
            Value to update by

        :param func:
            Optional arbitrary function that returns a float.  It should have
            the signature::

                func(existing_value, value)

            By default, func=np.add, so the cells will be simply incremented by
            `value`.

        :param cells:
            If True, assume that the distances `d1` and `d2` have already been
            normalized.
        """
        if not cells:
            d1 = self.normalize(d1)
            d2 = self.normalize(d2)
        super(HilbertNormalized, self).update(d1, d2, value, func)


class HilbertMatrix(HilbertNormalized):
    def __init__(self, file, genome, chrom, matrix_dim, incr_column=None,
                 default_chroms=True):
        """
        Subclass of HilbertNormalized that represents a genomic HilbertMatrix.

        If `default_chroms` is True, then only use the pybedtools-defined
        "default" chromosomes.  For example, this will be only the autosomes
        and X and Y for human, or just the euchromatic chromosomes for dm3.
        """
        self.file = file
        self.genome = genome
        self.chrom = chrom
        self.use_chrom_range = False

        # grab the dict of chrom lengths for this genome
        if isinstance(self.genome, basestring):
            self.chromdict = pbt.chromsizes(self.genome)
            if default_chroms:
                try:
                    self.chromdict = self.chromdict.default
                except AttributeError:
                    raise ValueError(
                        "Early version of pybedtools, or no chromosome "
                        "default set for genome %s.  Use "
                        "`default_chroms=False` instead." % self.genome)
        elif isinstance(self.genome, dict):
            self.chromdict = self.genome
        else:
            raise ValueError(
                '`genome` must be either a string assembly name '
                ' or a dictionary of chrom:(start, stop)')

        if self.chrom != "genome":
            chrom_range_tuple = get_interval_from_string(self.chrom)
            if chrom_range_tuple is None:
                # grab the length of the requested chromosome
                self.chrom_length = self.chromdict[self.chrom][1]
            else:
                (self.chrom,
                 self.range_start,
                 self.range_end) = chrom_range_tuple
                self.chrom_length = self.range_end - self.range_start
                self.use_chrom_range = True

            print self.chrom, "size: ",
        else:
            # using the entire genome for our coordinate system
            self.chrom_length = 0
            curr_offset = 0
            self.chrom_offsets = {}
            self.chrom_offsets_list = []
            self.chrom_names_list = []
            self.chrom_d = {}
            for chrom in self.chromdict:
                self.chrom_offsets[chrom] = curr_offset
                self.chrom_d[chrom] = curr_offset / (matrix_dim * matrix_dim)
                self.chrom_offsets_list.append(curr_offset)
                self.chrom_names_list.append(chrom)
                self.chrom_length += self.chromdict[chrom][1]
                curr_offset += self.chromdict[chrom][1]
            print "genome size: ",
        print self.chrom_length

        super(HilbertMatrix, self).__init__(matrix_dim, self.chrom_length)

        print "using matrix of size", self.matrix_dim, "there are", \
              self.ncells, "cells in the matrix and each cell represents", \
              int(self.dist_per_cell), "base pairs."

        self.incr_column = incr_column
        self.num_intervals = 0
        self.total_interval_length = 0
        chrom_offsets = []
        chrom_names = []
        self.temp_files = []

        # populate the matrix with the data contained in self.file
        self.build()
        self.dump_matrix()

    def _cleanup(self):
        for temp_file in self.temp_files:
            os.remove(temp_file)

    def _get_intervals(self):
        if not self.file.endswith('.bam'):
            main_file = pbt.BedTool(self.file)

            if not self.use_chrom_range:
                # no range requested, so just use the full file
                return main_file
            else:
                # range requested so intersect the main file
                # against the a range_file built from the range
                range_str = self.chrom + '\t' + \
                    str(self.range_start) + '\t' + \
                    str(self.range_end)
                range_file = pbt.BedTool(range_str, from_string=True)
                return main_file.intersect(range_file)
        else:
            # if we have a BAM file, we need to convert to
            # BEDGRAPH and force the use of the SCORE column
            # for incrementing the matrix
            bedg_filename = "." + self.file + ".bedg"
            bedg = open(bedg_filename, 'w')

            # use bedtools genomecov to create a BEDGRAPH
            # of the BAM file's coverage
            if self.chrom != "genome":
                cmd = "samtools view -b %s %s | \
                              bedtools genomecov -ibam - -bg" \
                              % (self.file, self.chrom)
            else:
                cmd = "bedtools genomecov -ibam %s -bg" % (self.file)

            # save the BEDGRAPH to a file for the Hilbert Curve
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            while True:
                line = proc.stdout.readline()
                if line != '':
                    bedg.write(line)
                else:
                    break
            bedg.close()

            # it's BEDGRAPH, so use the 4th column for the matrix cells
            self.incr_column = 4
            self.temp_files.append(bedg_filename)
            return pbt.BedTool(bedg_filename)

    def rc2chrom(self, row, col):
        """
        (chrom, start, stop) of a (row, col) coord
        """
        return self.get_chrom_range(col, row)

    def xy2chrom(self, x, y):
        """
        (chrom, start, stop) of an (x, y) coord
        """
        return self.get_chrom_range(x, y)

    def get_chrom_range(self, x, y):
        """
        Given an x,y coordinate in the matrix, compute the
        chrom, start and end coordinate that the cell represents.

        When plotting a single chromosome, this is merely a matter of
        getting the curve distance that the cell is from the origin
        and then multiplying that distance by the size (in bp) of
        each cell.

        For whole genome plots, we need to use the distance to figure out
        which chrom we are in and then from that, figure out how far we are
        into the chrom.
        """
        if self.chrom != "genome":
            matrix_dist = xy2d(self.matrix_dim, x, y)
            chrom = self.chrom
            start = matrix_dist * self.dist_per_cell
            end = start + self.dist_per_cell
        else:
            matrix_dist = xy2d(self.matrix_dim, x, y)
            bp_dist = matrix_dist * self.dist_per_cell

            idx = bisect.bisect_left(self.chrom_offsets_list, bp_dist) - 1
            chrom = self.chrom_names_list[idx]
            chrom_offset = self.chrom_offsets_list[idx]
            bp_dist_from_chrom_start = bp_dist - chrom_offset
            start = int(bp_dist_from_chrom_start / self.dist_per_cell) * \
                self.dist_per_cell
            end = start + self.dist_per_cell

        return chrom, int(start), int(end)

    def build(self):
        """
        Build the matrix.

        Need to compute a distance from the origin of the hilbert matrix start
        (0,0).  To do so, we create a normalization factor which is the length
        of the chrom divided by the number of cells in the matrix (d^2). Then,
        each coordinate is normalized by this constant to compute it's
        distance.  This distance is then converted to an x,y coordinate in the
        matrix using d2xy().
        """
        ivls = self._get_intervals()
        for ivl in ivls:
            self.num_intervals += 1
            self.total_interval_length += ivl.end - ivl.start

            start = ivl.start
            end = ivl.end

            # BED coords, which pybedtools uses under the hood, are zero-based
            # but do not include the stop position.
            #
            # The semantics of self.update() *do* include the stop position. So
            # subtract 1 from the end here -- but set it to the start if it's
            # a 1-bp feature.
            end = min(ivl.end - 1, ivl.start)

            if not self.chrom == "genome":
                if ivl.chrom != self.chrom:
                    continue
            else:
                offset = self.chrom_offsets[ivl.chrom]
                start = ivl.start + offset
                end = ivl.end + offset
            if self.incr_column is None:
                value = 1
            else:
                value = float(ivl[self.incr_column - 1])
            self.update(start, end, value=value)
        self._cleanup()

    def dump_matrix(self):
        """
        Export matrix as a text file in the same dir as the original
        """
        mat_dump = open(self.file + ".mtx", 'w')
        # header indicates the dimension of the matrix
        mat_dump.write(str(self.matrix_dim) + '\n')
        for r in xrange(self.matrix_dim):
            for c in xrange(self.matrix_dim):
                (chrom, start, end) = self.get_chrom_range(r, c)
                mat_dump.write('\t'.join(str(s) for s in [r, c,
                                         self.matrix[r][c],
                                         chrom, start,
                                         end]) + '\n')
        mat_dump.close()

    def norm_by_total_intervals(self):
        """
        Normalize matrix to total number of intervals
        """
        self.matrix /= self.num_intervals


class HilbertMatrixBigWig(HilbertMatrix):
    # Need to override build(), but otherwise just like a HilbertMatrix
    def __init__(self, *args, **kwargs):
        """
        Subclass of HilbertMatrix specifically for bigWig format files
        """
        super(HilbertMatrixBigWig, self).__init__(*args, **kwargs)

    def build(self):
        """
        Build the matrix.

        Since bigWig files are essentially pre-summarized, this just extracts
        the chrom/start/stop represented by each cell in the matrix and fills
        it with the value from the bigWig file.
        """
        self.bigwig = BigWigFile(open(self.file))

        chrom_rc, chrom_bins = self.chrom2rc()

        if self.chrom == 'genome':
            chroms = self.chromdict.keys()

        else:
            chroms = [self.chrom]

        for chrom in chroms:
            rc = chrom_rc[chrom]
            nbins = chrom_bins[chrom]

            start, stop = self.chromdict[chrom]
            results = self.bigwig.summarize(chrom, start, stop, nbins)
            values = results.sum_data / results.valid_count
            values[np.isnan(values)] = 0

            self.matrix[rc[:,0], rc[:, 1]] = values

        self._cleanup()


    def chrom2rc(self):
        """
        Return a dictionary of {chrom: (rows, cols)} and {chrom: nbins}
        """
        precomputed = np.load(
            os.path.join(
                os.path.dirname(__file__),
                'precomputed.npz'))
        rc = precomputed['_%s' % self.matrix_dim]

        d = {}
        bins = {}
        last_stop = 0
        for chrom, startstop in self.chromdict.items():
            start, stop = startstop
            frac = self.chromdict[chrom][1] / float(self.chrom_length)
            nbins = int(frac * (self.matrix_dim * self.matrix_dim))
            d_start = last_stop
            d_stop = d_start + nbins
            d[chrom] = rc[d_start:d_stop, :]
            bins[chrom] = nbins
            last_stop += nbins
        return d, bins