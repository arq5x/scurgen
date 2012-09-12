import numpy as np
import pybedtools
from scurgen import hilbert
from scurgen.plotting import debug_plot
from nose.tools import assert_raises


def test_power_of_2():
    h = hilbert.HilbertBase(16)
    assert_raises(ValueError, hilbert.HilbertBase, 15)


def test_size():
    h = hilbert.HilbertBase(16)
    assert(h.matrix.shape == (16, 16))
    assert(h.matrix.sum() == 0)


def test_invariant_normalize():
    h = hilbert.HilbertNormalized(16, length=1000)
    # had to change:
    #   self.norm_factor = self.length / self.ncells
    #
    # to:
    #   self.norm_factor = self.length / float(self.ncells)
    #
    # to get this to pass
    assert h.normalize(h.length) == h.ncells

    # Test support for small lengths
    h = hilbert.HilbertNormalized(16, length=0.5)
    assert h.normalize(h.length) == h.ncells


def test_normalize():
    h = hilbert.HilbertNormalized(16, length=16 * 16)
    assert h.norm_factor == 1
    assert h.normalize(1) == 1

    # This one should behave identically to HilbertBase of the same size
    hb = hilbert.HilbertBase(16)
    h.update(0, 3)
    hb.update(0, 3)
    assert np.all(hb.matrix == h.matrix)

    # cells are now half the size
    h = hilbert.HilbertNormalized(16, length=0.5 * 16 * 16)
    assert h.norm_factor == 0.5
    assert h.normalize(0.5) == 1
    assert h.normalize(1) == 2


def test_cell_fill():
    """
    makes sure the right (rows, cols) get filled as expected when using cells
    as distance
    """
    h = hilbert.HilbertBase(16)

    # append 1 to the first cell
    h.update(0, 0)
    assert h.matrix[0][0] == 1
    assert h.matrix.sum() == 1

    # append 5 to the first 2 cells
    h.update(0, 1, value=5)

    # can visually check what it should look like
    #debug_plot(h)

    assert h.matrix.sum() == 11
    assert h.matrix[0][0] == 6
    assert h.matrix[1][0] == 5

    # arbitrary function
    def overwrite(orig, new):
        if orig > 0:
            return 0
        return new

    # apply arbitrary function
    h.update(0, 2, value=99, func=overwrite)

    #debug_plot(h)
    assert h.matrix[1][1] == 99

    # everything from before should have been set back to zero
    assert h.matrix.sum() == 99


def test_norm_dist_fill():
    """
    make sure the right (rows, cols) get filled when using distance rather than
    cell count
    """
    # test a variety of lengths
    lengths = [0.5, 1, 16, 1000]
    cells_to_fill = 3
    for length in lengths:
        h = hilbert.HilbertNormalized(16, length)

        # minus 1 because in the h.update() below, the start, 0, counts as the
        # first
        dist_to_fill = (cells_to_fill - 1) * h.dist_per_cell

        h.update(0, dist_to_fill)

        assert h.matrix[0, 0] == 1
        assert h.matrix[1, 0] == 1
        assert h.matrix[1, 1] == 1
        assert h.matrix.sum() == 3

        h.reset()

        h.update(1 * h.dist_per_cell, dist_to_fill)

        assert h.matrix[0, 0] == 0
        assert h.matrix[1, 0] == 1
        assert h.matrix[1, 1] == 1
        assert h.matrix.sum() == 2


def test_xy2d_d2xy():
    n = 8
    for d in range(n * n):
        x, y = hilbert.d2xy(n, d)
        d2 = hilbert.xy2d(n, x, y)
        assert d == d2

def test_masked():
    h = hilbert.HilbertBase(16)
    h.update(0, 5, value=10)
    h.update(0, 1, value=50)

    assert h.matrix.sum() == 160
    assert h.masked.sum() == 160

    h.mask_low_values(min_val=10)

    # masking should not change the original data
    assert h.matrix.sum() == 160

    # but the masked one has changed
    assert h.masked.sum() == 120

    # mask should have only masked 2 cells
    assert (h.masked.mask == False).sum() == 2



def test_chrom_range():
    dim = 8

    # Create a BedTool object that has features every 1 bp through the end of
    # the chrom, and scores continuously increment.
    def generator():
        for i in range(dim * dim):
            yield pybedtools.create_interval_from_list([
                'chr1', str(i), str(i+1), str(i)])

    one_bp_bt = pybedtools.BedTool(generator()).saveas()

    # conveniently-sized genome completely encompasses intervals in one_bp_bt
    genome = {'chr1': (0, dim * dim)}

    h = hilbert.HilbertMatrix(
            # provide just the filename of the BedTool
            file=one_bp_bt.fn,
            genome=genome,
            chrom='chr1',
            matrix_dim=dim,
            # incr_column is 1-based!
            incr_column=4)

    # Everything in the matrix should contain all the scores in the bedgraph
    # file, and given how everything was constructed it ought to be the same as
    # the sum of integers from 0 to dim^2
    assert h.matrix.sum() == sum([int(i[-1]) for i in one_bp_bt]) == sum(range(dim * dim))

    #debug_plot(h)

    # be paranoid to make sure we're using the right matrix...this checks that
    # the values in the cells are correct...
    assert h.matrix[0, 0] == 0
    assert h.matrix[0, 1] == 1
    assert h.matrix[1, 1] == 2
    assert h.matrix[1, 0] == 3
    assert h.matrix[2, 0] == 4
    assert h.matrix[3, 0] == 5

    # ...and this checks that those same cells return the correct chrom_range
    assert h.get_chrom_range(0, 0) == ('chr1', 0, 1)
    assert h.get_chrom_range(0, 1) == ('chr1', 1, 2)
    assert h.get_chrom_range(1, 1) == ('chr1', 2, 3)
    assert h.get_chrom_range(1, 0) == ('chr1', 3, 4)
    assert h.get_chrom_range(2, 0) == ('chr1', 4, 5)
    assert h.get_chrom_range(3, 0) == ('chr1', 5, 6)
