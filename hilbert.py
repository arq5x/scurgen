import numpy as np
import pybedtools as pbt
import matplotlib
from pybedtools import genome_registry

def rot(n, x, y, rx, ry):
    if (ry == 0):
        if (rx == 1):
            x = n-1 - x
            y = n-1 - y
        # Swap x and y
        x, y = y, x

def d2xy(n, d):
    t = d
    x = y = 0
    
    s = 1
    while s < n:
        rx = 1 & (t/2)
        ry = 1 & (t ^ rx)
        rot(s, x, y, rx, ry)
        x += (s * rx)
        y += (s * ry)
        t = t / 4
        s *= 2
    return (x,y)


class HilbertMatrix(object):
    
    def __init__(self, file, genome, chrom, matrix_dim, incr_column=-1):
        self.file = file
        self.genome = genome
        self.chrom = chrom
        
        # grab the dict of chrom lengths for this genome
        self.chromdict = pbt.chromsizes(self.genome)
        # grab the length of the requested genome
        self.chrom_length = self.chromdict[self.chrom][1]
        self.m_dim = matrix_dim
        self.cells = self.m_dim * self.m_dim
        self.norm_factor = int(self.chrom_length / self.cells)
        self.incr_column = incr_column
        
        chromdict = pbt.chromsizes(genome)
        # populate the matrix with the data contained in self.file
        self.build()

    def _update_matrix(self, coords, increment=1):
        x = coords[0]
        y = coords[1]
        self.matrix[x][y] += increment
    
    def build(self):
    
        # must compute a distance from the 
        # origing of the hilbert matrix start (0,0.
        # to do so, we create a normalization factor which
        # is the length of the chrom divided by the number
        # of cells in the matrix (d^2). then, each coordinate
        # "normalized by this constant to compute it's distance.
        # this distance is then converted to an x,y coordinate
        # in the matrix using d2xy
        
        # initialize the matrix
        self.matrix = np.zeros((self.m_dim,self.m_dim), dtype=np.float)
        for ivl in pbt.BedTool(self.file):
        
            start_dist = int(ivl.start / self.norm_factor)
            end_dist   = int(ivl.end / self.norm_factor)

            for dist in xrange(start_dist, end_dist + 1):
                coords = d2xy(self.m_dim, dist)
                if self.incr_column:
                    self._update_matrix(coords)
                else:
                    self._update_matrix(coords, \
                        increment=float(ivl[self.incr_column]))
                        
    def mask_zeros(self):
        rows, cols = self.matrix.shape
        for r in range(rows):
            for c in range(cols):
                if self.matrix[r][c] == 0:
                    self.matrix[r][c] = np.NaN
