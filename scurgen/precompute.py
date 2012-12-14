"""
Pre-compute distance to row, col for a range of powers-of-2.
"""
from scurgen.hilbert import d2rc
import numpy as np
import sys
d = {}
for order in range(2, 11):
    dim = 2 ** order
    n = dim ** 2
    print 'computing order=%s (dim=%s)' % (order, dim)
    sys.stdout.flush()
    d['_%s' % dim] = np.array([d2rc(n, i) for i in range(n)])

np.savez_compressed('precomputed', **d)
