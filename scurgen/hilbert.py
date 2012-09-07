import numpy as np
import pybedtools as pbt
import sys
import os
import matplotlib
import subprocess
from pybedtools import genome_registry


def rot(n, x, y, rx, ry):
    if (ry == 0):
        if (rx == 1):
            x = n - 1 - x
            y = n - 1 - y
        # Swap x and y
        x, y = y, x
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
    return d


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


class HilbertBase(object):
    def __init__(self, m_dim):
        """
        Initializes a basic Hilbert curve.

        :param m_dim:
            The number of dimensions on a side (power of 2)
        """
        if np.log2(m_dim) % 1 != 0:
            raise ValueError('m_dim %s not a power of 2' % m_dim)
        self.m_dim = m_dim
        self.ncells = m_dim * m_dim
        self.matrix = np.zeros((self.m_dim, self.m_dim), dtype=np.float)

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
            x, y = d2xy(self.m_dim, dist)
            self.matrix[x, y] = func(self.matrix[x, y], value)

    def curve(self):
        """
        Returns a 3-tuple of the x-coords, y-coords, and labels for the curve.
        The curve is rotated such that (0,0) is in the upper left and (-1,-1)
        is in the lower left in order to be consistent with matplotlib's
        imshow() default axes origin.
        """
        xs, ys = [], []
        for i in range(self.ncells):
            x, y = d2xy(self.m_dim, i)
            xi = y
            yi = self.m_dim - x - 1
            xs.append(xi)
            ys.append(yi)
        return np.array(xs), np.array(ys), range(self.ncells)[::-1]


class HilbertNormalized(HilbertBase):
    def __init__(self, m_dim, length):
        """
        Hilbert curve class that handles distance in arbitrary units (e.g.,
        genomic bp) instead of in cell numbers.

        :param m_dim:
            The number of dimensions on a side (power of 2)

        :param length:
            The total length represented by this curve
        """
        super(HilbertNormalized, self).__init__(m_dim)
        self.length = length
        if self.ncells >= length:
            raise ValueError("too many cells for this length")
        self.norm_factor = self.length / self.ncells

    def normalize(self, d):
        """
        Convert distance from bp to cell number
        """
        return int(d / self.norm_factor)

    def update(self, d1, d2, value=1, func=np.add, cells=False):
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
        else:
            # using the entire genome for our coordinate system
            self.chrom_length = 0
            curr_offset = 0
            self.chrom_offsets = {}
            for chrom in self.chromdict:
                self.chrom_offsets[chrom] = curr_offset
                self.chrom_length += self.chromdict[chrom][1]
                curr_offset += self.chromdict[chrom][1]

        self.incr_column = incr_column
        self.num_intervals = 0
        self.total_interval_length = 0
        chromdict = pbt.chromsizes(genome)
        self.temp_files = []

        super(HilbertMatrix, self).__init__(matrix_dim, self.chrom_length)

        # populate the matrix with the data contained in self.file
        self.build()

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

    def mask_low_values(self, min_val=0):
        rows, cols = self.matrix.shape
        for r in range(rows):
            for c in range(cols):
                if self.matrix[r][c] <= min_val:
                    self.matrix[r][c] = np.NaN

    def norm_by_total_intervals(self):
        rows, cols = self.matrix.shape
        for r in range(rows):
            for c in range(cols):
                self.matrix[r][c] /= self.num_intervals
