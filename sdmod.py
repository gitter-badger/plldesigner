"""
Implementation of a sigma delta modulator and the functions related
"""

import numpy as np
from numpy import (zeros, arange, log10, sin, pi)
import matplotlib.pyplot as plt


class SDModulator(object):
    """ Creates a Sigma-Delta modulator object

    Parameters
    ==========
    modtype : string
        type of SDModulator implemented types are:
        ('mash')
    *arg  : Positional arguments for the different types of SDMOD
    **argk :
    """
    def __init__(self, modtype, *arg, **argk):
        """ """
        self.modtype = modtype
        func = {'mash': gen_mash}
        self.seq, self.cycles = func[modtype](*arg, **argk)

    def plotseq(self, *arg, **argk):
        xval = np.arange(len(self.seq))
        plt.step(xval, self.seq, *arg, **argk)

    def LdB_theoretical(self, fm, fref, N=1.0):
        func = {'mash': L_mash_dB}
        return func[self.modtype](fm, fref, N)


def gen_mash(order, N, K, init=()):
    """ Generates a mash type $\Sigma-\Delta$ sequence

    Parameters
    ----------
    order : int
        order of the $\Sigma-\Delta$ modulator.
    N : int
        Number of bits of the modulator.
    K : int
        Value being converted
    init: tuple
        init is a tuple that initialize the register of the
        mash sd, at with length equal to the order.

    Returns
    -------
    sd :
    cycles : int
        Period of the sequence

    Commentaries
    ------------
    This implementation is not really effective from the computational
    point of view but is representative of the architecture of the
    system.
    """

    # Assertions
    assert len(init) == order or len(init) == 0, (
        'Initial condition length must be equal to the order of the modulator')
    assert order >= 1 and order <= 4, 'This orders have not been implemented'

    # Modulator of order 1
    _maxvalue = 2**N-1
    L = len(K)
    if order == 1:
        # initialize the registers
        over0 = zeros(L, dtype=np.int)
        over0[0] = 1
        stat0 = zeros(L, dtype=np.int)
        if len(init) == 1:
            stat0[0] = init[0]
        for j in arange(1, L):
            stat0[j] = stat0[j-1] + K[j-1]
            if stat0[j] > _maxvalue:
                over0[j] = 1
                stat0[j] -= _maxvalue + 1
        sd = over0
        cycles = np.where(stat0 == 0)[0]
        if len(cycles) > 1:
            cycles = cycles[1]-cycles[0]
        else:
            cycles = -1

    # Modulator of order 2
    elif order == 2:
        # initialize the registers
        stat0 = zeros(L, dtype=np.int)
        stat1 = zeros(L, dtype=np.int)
        over0 = zeros(L, dtype=np.int)
        over1 = zeros(L, dtype=np.int)
        if len(init) == 2:
            stat0[0] = init[0]
            stat1[0] = init[1]

        # Implement the SDM
        for j in arange(1, L):
            stat0[j] = stat0[j-1] + K[j-1]
            if stat0[j] > _maxvalue:
                over0[j] = 1
                stat0[j] -= _maxvalue + 1
            stat1[j] = stat1[j-1] + stat0[j]
            if stat1[j] > _maxvalue:
                over1[j] = 1
                stat1[j] -= _maxvalue + 1
        sd = over0 + over1 - np.hstack(([0], over1[:-1]))
        stat = stat0 + stat1
        cycles = np.where(stat == 0)[0]
        if len(cycles) > 1:
            cycles = cycles[1] - cycles[0]
        else:
            cycles = -1

    # Modulator of order 3
    elif order == 3:
        # initialize the registers
        stat0 = zeros(L, dtype=np.int)
        stat1 = zeros(L, dtype=np.int)
        stat2 = zeros(L, dtype=np.int)
        if len(init) == 3:
            stat0[0] = init[0]
            stat1[0] = init[1]
            stat2[0] = init[2]
        over0 = zeros(L, dtype=np.int)
        over0[0] = 1
        over1 = zeros(L, dtype=np.int)
        over2 = zeros(L, dtype=np.int)

        # Implement the SDM
        for j in arange(1, L):
            stat0[j] = stat0[j-1] + K[j-1]
            if stat0[j] > _maxvalue:
                over0[j] = 1
                stat0[j] -= _maxvalue + 1

            stat1[j] = stat1[j-1] + stat0[j]
            if stat1[j] > _maxvalue:
                over1[j] = 1
                stat1[j] -= _maxvalue + 1

            stat2[j] = stat2[j-1] + stat1[j]
            if stat2[j] > _maxvalue:
                over2[j] = 1
                stat2[j] -= _maxvalue + 1
        sd = over0
        sd += over1 - np.hstack(([0], over1[:-1]))
        sd += over2 - 2*np.hstack(([0], over2[:-1]))
        sd += np.hstack(([0], [0], over2[:-2]))

        stat = stat0 + stat1 + stat2
        cycles = np.where(stat == 0)[0]
        if len(cycles) > 1:
            cycles = cycles[1] - cycles[0]
        else:
            cycles = -1

    elif order == 4:
        # initialize the registers
        stat0 = zeros(L, dtype=np.int)
        stat1 = zeros(L, dtype=np.int)
        stat2 = zeros(L, dtype=np.int)
        stat3 = zeros(L, dtype=np.int)
        if len(init) == 4:
            stat0[0] = init[0]
            stat1[0] = init[1]
            stat2[0] = init[2]
            stat3[0] = init[4]
        over0 = zeros(L, dtype=np.int)
        over0[0] = 1
        over1 = zeros(L, dtype=np.int)
        over2 = zeros(L, dtype=np.int)
        over3 = zeros(L, dtype=np.int)
        # Implement the SDM
        for j in arange(1, L):
            stat0[j] = stat0[j-1] + K[j-1]
            if stat0[j] > _maxvalue:
                over0[j] = 1
                stat0[j] -= _maxvalue + 1

            stat1[j] = stat1[j-1] + stat0[j]
            if stat1[j] > _maxvalue:
                over1[j] = 1
                stat1[j] -= _maxvalue + 1

            stat2[j] = stat2[j-1] + stat1[j]
            if stat2[j] > _maxvalue:
                over2[j] = 1
                stat2[j] -= _maxvalue + 1
            stat3[j] = stat3[j-1] + stat2[j]
            if stat3[j] > _maxvalue:
                over3[j] = 1
                stat3[j] -= _maxvalue + 1

        sd = (over0 +
               over1 - np.hstack(([0], over1[:-1])) +
               over2 - 2*np.hstack(([0], over2[:-1])) +
               np.hstack(([0], [0], over2[:-2])) +
               over3 - 3*np.hstack(([0], over3[:-1])) +
               3*np.hstack(([0], [0], over3[:-2])) -
               np.hstack(([0], [0], [0], over3[:-3]))
        )

        stat = stat0 + stat1 + stat2 + stat3
        cycles = np.where(stat == 0)[0]
        if len(cycles) > 1:
            cycles = cycles[1] - cycles[0]
        else:
            cycles = -1

    return sd, cycles

def L_mash_dB(m, fm, fref, N=1.0):
    """ Phase noise theoretical value of noise produced by a mash111 SDM

    This procedure calculates the noise at the output of the SD modulator

    Parameters
    ----------
    m : int
        The order of the SD modulator
    fm : array_like
        Frequency offsets were the  noise is calculated
    fref : float
        Reference frequency that is used to compare the output of the SD
        modulator
    N: float
        It is the average sdision ratio.
    """
    return 10*log10((2*pi)**2/(12*fref)*(2*sin(pi*fm/fref))**(2*(m-1))/N**2)

if __name__ == "__main__":
    from numpy.testing import assert_almost_equal
    import numpy.random as rnd
    # Test order one assert the mean value
    floatnum = rnd.rand()*np.ones(100000)

    # order one
    sequence, cycles = gen_mash(1, 19, (floatnum*2**19).astype(int))
    assert_almost_equal(sequence.mean(), floatnum.mean(), 4)
    # order two
    sequence, cycles = gen_mash(2, 19, (floatnum*2**19).astype(int))
    assert_almost_equal(sequence.mean(), floatnum.mean(), 4)

    # order three
    sequence, cycles = gen_mash(3, 19, (floatnum*2**19).astype(int))
    assert_almost_equal(sequence.mean(), floatnum.mean(), 4)

    # order three
    sequence, cycles = gen_mash(3, 19, 0.25*np.ones(100000)*2**19)
    assert_almost_equal(sequence.mean(), 0.25, 4)

    # order three
    sequence, cycles = gen_mash(4, 19, 0.25*np.ones(100000)*2**19)
    assert_almost_equal(sequence.mean(), 0.25, 4)

    # Using the class the test is done as:
    sd_mash = SDModulator('mash', 3, 19,
                    (np.array([0.323232]*100000)*2**19).astype(int))
    assert_almost_equal(sd_mash.seq.mean(), 0.323232, 4)
