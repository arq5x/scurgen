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
            x = n-1 - x
            y = n-1 - y
        # Swap x and y
        x, y = y, x
    return (x, y)
        
def d2xy(n, d):
    t = d
    x = y = 0
    
    s = 1
    while s < n:
        rx = 1 & (t/2)
        ry = 1 & (t ^ rx)
        (x, y) = rot(s, x, y, rx, ry)
        x += (s * rx)
        y += (s * ry)
        t = t / 4
        s *= 2
    return (x,y)


class HilbertMatrix(object):
    
    def __init__(self, file, genome, chrom, matrix_dim, incr_column=None):
        self.file = file
        self.genome = genome
        self.chrom = chrom
        
        # grab the dict of chrom lengths for this genome
        self.chromdict = pbt.chromsizes(self.genome)
        
        if self.chrom != "genome":
            # grab the length of the requested genome
            self.chrom_length = self.chromdict[self.chrom][1]
            print self.chrom, "size: ", 
        else:
            # using the entire genome for our coordinate system
            self.chrom_length = 0
            curr_offset = 0
            self.chrom_offsets = {}
            for chrom in self.chromdict:
                self.chrom_offsets[chrom] = curr_offset
                self.chrom_length += self.chromdict[chrom][1]
                curr_offset += self.chromdict[chrom][1]
            print "genome size: ",
        print self.chrom_length
        
        self.m_dim = matrix_dim
        self.cells = self.m_dim * self.m_dim
        self.norm_factor = int(self.chrom_length / self.cells)
        
        print "using matrix of size", self.m_dim, "there are", \
              self.cells, "cells in the matrix and each cell represents", \
              self.norm_factor, "base pairs."
        
        self.incr_column = incr_column
        self.num_intervals = 0
        self.total_interval_length = 0
        chromdict = pbt.chromsizes(genome)
        self.temp_files = []

        # populate the matrix with the data contained in self.file
        self.build()

    def _cleanup(self):
        for temp_file in self.temp_files:
            os.remove(temp_file)

    def _update_matrix(self, coords, increment=1):
        x = coords[0]
        y = coords[1]
        self.matrix[x][y] += increment
        
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
        # initialize the matrix
        self.matrix = np.zeros((self.m_dim,self.m_dim), dtype=np.float)
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
                end   = ivl.end + offset
            
                
            # figure out what cell the start and end coords
            # of the interval belong in.
            # most of the time, the interval will fit in a single cell
            start_dist = int(start / self.norm_factor)
            end_dist   = int(end / self.norm_factor)
            
            # however, we must populate EVERY cell that the 
            # interval spans.
            for dist in xrange(start_dist, end_dist + 1):
                coords = d2xy(self.m_dim, dist)
                if self.incr_column is None:
                    self._update_matrix(coords)
                else:
                    self._update_matrix(coords, \
                        increment=float(ivl[self.incr_column - 1]))
        self._cleanup()
        
    def mask_low_values(self, min_val = 0):
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
