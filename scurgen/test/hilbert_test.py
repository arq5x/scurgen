from scurgen import hilbert
from nose.tools import assert_raises

def test_power_of_2():
    h = hilbert.HilbertBase(16)
    assert_raises(ValueError, hilbert.HilbertBase, 15)

def test_size():
    h = hilbert.HilbertBase(16)
    assert(h.matrix.shape == (16, 16))
    assert(h.matrix.sum() == 0)
