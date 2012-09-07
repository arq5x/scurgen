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


def test_cell_fill():
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
