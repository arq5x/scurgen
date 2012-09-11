import numpy as np
import pybedtools as pbt
import sys
import os
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
    x = y = 0

    s = 1
    while s < n:
        rx = 1 & (t / 2)
        ry = 1 & (t ^ rx)
        (x, y) = rot(s, x, y, rx, ry)
        x += (s * rx)
        y += (s * ry)
        t = t / 4
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
        rot(s, x, y, rx, ry)
        s /= 2
    return d


class Interval(object):

    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end


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
            (self.matrix_dim, self.matrix_dim), dtype=np.float)

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
            x, y = d2xy(self.matrix_dim, dist)
            self.matrix[x, y] = func(self.matrix[x, y], value)

    def reset(self):
        """
        Resets the matrix to zeros everywhere
        """
        self.matrix[:] = 0.0

    def curve(self):
        """
        Returns a 3-tuple of the x-coords, y-coords, and labels for the curve.
        The (x, y) coords correspond to (col, row) positions in a matrix such
        that that the following will work::

            x, y, labels = x.curve()
            plt.imshow(x.matrix)
            plt.plot(x, y)

        In this case the (x, y) coords fall within the center of each cell.
        """
        # Note: there's probably some more elegant way to do this but this
        # seems to work...
        xs, ys = [], []
        for i in range(self.ncells):
            x, y = d2xy(self.matrix_dim, i)
            xi = y
            yi = self.matrix_dim - x - 1
            xs.append(xi)
            ys.append(yi)

        # Reverse everything
        return np.array(xs)[::-1], np.array(ys)[::-1], range(self.ncells)


class HilbertNormalized(HilbertBase):
    def __init__(self, matrix_dim, length):
        """
        Hilbert curve class that handles distance in arbitrary units (e.g.,
        genomic bp) instead of in cell numbers.

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
    def __init__(self, file, genome, chrom, matrix_dim, incr_column=None):
        self.file = file
        self.genome = genome
        self.chrom = chrom

        # grab the dict of chrom lengths for this genome
        self.chromdict = pbt.chromsizes(self.genome)

        if self.chrom != "genome":
            # grab the length of the requested chromosome
            self.chrom_length = self.chromdict[self.chrom][1]
            print self.chrom, "size: ",
        else:
            # using the entire genome for our coordinate system
            self.chrom_length = 0
            curr_offset = 0
            self.chrom_offsets = {}
            self.chrom_offsets_list = []
            self.chrom_names_list = []
            for chrom in self.chromdict:
                self.chrom_offsets[chrom] = curr_offset
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
        chromdict = pbt.chromsizes(genome)
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
            return pbt.BedTool(self.file)
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

    def get_chrom_range(self, x, y):
        """
        Given an x,y coordinate in the matrix, compute the
        chrom, start and end coordinate tht the cell represents.

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
        must compute a distance from the
        origing of the hilbert matrix start (0,0.
        to do so, we create a normalization factor which
        is the length of the chrom divided by the number
        of cells in the matrix (d^2). then, each coordinate
        "normalized by this constant to compute it's distance.
        this distance is then converted to an x,y coordinate
        in the matrix using d2xy
        """
        ivls = self._get_intervals()
        for ivl in ivls:
            self.num_intervals += 1
            self.total_interval_length += ivl.end - ivl.start

            start = ivl.start
            end = ivl.end
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

    def mask_low_values(self, min_val=0):
        self.matrix = np.ma.masked_where(self.matrix <= min_val, self.matrix)

    def norm_by_total_intervals(self):
        rows, cols = self.matrix.shape
        for r in range(rows):
            for c in range(cols):
                self.matrix[r][c] /= self.num_intervals
