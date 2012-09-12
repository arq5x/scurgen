import numpy as np
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
