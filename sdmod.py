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
        self.seq, self.per = func[modtype](*arg, **argk)

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
    div :
    per : int
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
    MAXVAL = 2**N-1
    L = len(K)
    if order == 1:
        # initialize the registers
        over1 = zeros(L, dtype=np.int)
        over1[0] = 1
        stat1 = zeros(L, dtype=np.int)
        if len(init) == 1:
            stat1[0] = init[0]
        for j in arange(1, L):
            stat1[j] = stat1[j-1] + K[j-1]
            if stat1[j] > MAXVAL:
                over1[j] = 1
                stat1[j] = stat1[j] - (MAXVAL+1)
        div = over1
        per = np.where(stat1 == 0)[0]
        if len(per) > 1:
            per = per[1]-per[0]
        else:
            per = -1

    # Modulator of order 2
    elif order == 2:
        # initialize the registers
        stat1 = zeros(L, dtype=np.int)
        stat2 = zeros(L, dtype=np.int)
        over1 = zeros(L, dtype=np.int)
        over2 = zeros(L, dtype=np.int)
        if len(init) == 2:
            stat1[0] = init[0]
            stat2[0] = init[1]

        # Implement the SDM
        for j in arange(1, L):
            stat1[j] = stat1[j-1] + K[j-1]
            if stat1[j] > MAXVAL:
                over1[j] = 1
                stat1[j] -= MAXVAL + 1
            stat2[j] = stat2[j-1] + stat1[j]
            if stat2[j] > MAXVAL:
                over2[j] = 1
                stat2[j] -= MAXVAL + 1
        div = over1 + over2 - np.hstack(([0], over2[:-1]))
        stat = stat1 + stat2
        per = np.where(stat == 0)[0]
        if len(per) > 1:
            per = per[1] - per[0]
        else:
            per = -1

    # Modulator of order 3
    elif order == 3:
        # initialize the registers
        stat1 = zeros(L, dtype=np.int)
        stat2 = zeros(L, dtype=np.int)
        stat3 = zeros(L, dtype=np.int)
        if len(init) == 3:
            stat1[0] = init[0]
            stat2[0] = init[1]
            stat3[0] = init[2]
        over1 = zeros(L, dtype=np.int)
        over1[0] = 1
        over2 = zeros(L, dtype=np.int)
        over3 = zeros(L, dtype=np.int)

        # Implement the SDM
        for j in arange(1, L):
            stat1[j] = stat1[j-1] + K[j-1]
            if stat1[j] > MAXVAL:
                over1[j] = 1
                stat1[j] -= MAXVAL + 1

            stat2[j] = stat2[j-1] + stat1[j]
            if stat2[j] > MAXVAL:
                over2[j] = 1
                stat2[j] -= MAXVAL + 1

            stat3[j] = stat3[j-1] + stat2[j]
            if stat3[j] > MAXVAL:
                over3[j] = 1
                stat3[j] -= MAXVAL + 1
        div = over1
        div += over2 - np.hstack(([0], over2[:-1]))
        div += over3 - 2*np.hstack(([0], over3[:-1]))
        div += np.hstack(([0], [0], over3[:-2]))

        stat = stat1 + stat2 + stat3
        per = np.where(stat == 0)[0]
        if len(per) > 1:
            per = per[1] - per[0]
        else:
            per = -1

    elif order == 4:
        # initialize the registers
        stat1 = zeros(L, dtype=np.int)
        stat2 = zeros(L, dtype=np.int)
        stat3 = zeros(L, dtype=np.int)
        stat4 = zeros(L, dtype=np.int)
        if len(init) == 4:
            stat1[0] = init[0]
            stat2[0] = init[1]
            stat3[0] = init[2]
            stat4[0] = init[4]
        over1 = zeros(L, dtype=np.int)
        over1[0] = 1
        over2 = zeros(L, dtype=np.int)
        over3 = zeros(L, dtype=np.int)
        over4 = zeros(L, dtype=np.int)
        # Implement the SDM
        for j in arange(1, L):
            stat1[j] = stat1[j-1] + K[j-1]
            if stat1[j] > MAXVAL:
                over1[j] = 1
                stat1[j] -= MAXVAL + 1

            stat2[j] = stat2[j-1] + stat1[j]
            if stat2[j] > MAXVAL:
                over2[j] = 1
                stat2[j] -= MAXVAL + 1

            stat3[j] = stat3[j-1] + stat2[j]
            if stat3[j] > MAXVAL:
                over3[j] = 1
                stat3[j] -= MAXVAL + 1
            stat4[j] = stat4[j-1] + stat3[j]
            if stat4[j] > MAXVAL:
                over4[j] = 1
                stat4[j] -= MAXVAL + 1

        div =  (over1 +
                over2 - np.hstack(([0], over2[:-1])) +
                over3 - 2*np.hstack(([0], over3[:-1])) +
                np.hstack(([0], [0], over3[:-2])) +
                over4 - 3*np.hstack(([0], over4[:-1])) +
                3*np.hstack(([0], [0], over4[:-2])) -
                np.hstack(([0], [0], [0], over4[:-3]))
                )


        stat = stat1 + stat2 + stat3 + stat4
        per = np.where(stat == 0)[0]
        if len(per) > 1:
            per = per[1] - per[0]
        else:
            per = -1

    return div, per


def L_mash_dB(m, fm, fref, N=1.0):
    """ Phase noise theoretical value of noise produced by a mash SDM

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
        It is the average division ratio.
    """
    return 10*log10((2*pi)**2/(12*fref)*(2*sin(pi*fm/fref))**(2*(m-1))/N**2)

if __name__ == "__main__":
    from numpy.testing import assert_almost_equal
    import numpy.random as rnd
    # Test order one assert the mean value
    floatnum = rnd.rand()*np.ones(100000)

    # order one
    sequence, per = gen_mash(1, 19, (floatnum*2**19).astype(int))
    assert_almost_equal(sequence.mean(), floatnum.mean(), 4)
    # order two
    sequence, per = gen_mash(2, 19, (floatnum*2**19).astype(int))
    assert_almost_equal(sequence.mean(), floatnum.mean(), 4)

    # order three
    sequence, per = gen_mash(3, 19, (floatnum*2**19).astype(int))
    assert_almost_equal(sequence.mean(), floatnum.mean(), 4)

    # order three
    sequence, per = gen_mash(3, 19, 0.25*np.ones(100000)*2**19)
    assert_almost_equal(sequence.mean(), 0.25, 4)

    # order three
    sequence, per = gen_mash(4, 19, 0.25*np.ones(100000)*2**19)
    assert_almost_equal(sequence.mean(), 0.25, 4)

    # Using the class the test is done as:
    sd_mash = SDModulator('mash', 3, 19,
                    (np.array([0.323232]*100000)*2**19).astype(int))
    assert_almost_equal(sd_mash.seq.mean(), 0.323232, 4)
