"""
Implementation of a sigma delta modulator and the functions related
"""

import numpy as np
from numpy import (zeros, arange, log10, sin, pi)
import matplotlib.pyplot as plt
from .pnoise import Pnoise


class SDModulator(object):
    """ Creates a Sigma-Delta modulator object

    Parameters
    ----------
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

    def LdB_theoretical(self, fm, fref, n=1.0):
        func = {'mash': L_mash_dB}
        return func[self.modtype](fm, fref, n)


def gen_mash(order, n, k, init=()):
    """ Generates a mash type $\Sigma-\Delta$ sequence

    Parameters
    ----------
    order : int
        order of the $\Sigma-\Delta$ modulator.
    n : int
        Number of bits of the modulator.
    k : int
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
    assert 1 <= order <= 4, 'This orders have not been implemented'

    # Modulator of order 1
    _maxvalue = 2 ** n - 1
    L = len(k)
    if order == 1:
        # initialize the registers
        overflow0 = zeros(L, dtype=np.int)
        overflow0[0] = 1
        state0 = zeros(L, dtype=np.int)
        if len(init) == 1:
            state0[0] = init[0]
        for j in arange(1, L):
            state0[j] = state0[j - 1] + k[j - 1]
            if state0[j] > _maxvalue:
                overflow0[j] = 1
                state0[j] -= _maxvalue + 1
        sd = overflow0
        cycles = np.where(state0 == 0)[0]
        if len(cycles) > 1:
            cycles = cycles[1] - cycles[0]
        else:
            cycles = -1

    # Modulator of order 2
    elif order == 2:
        # initialize the registers
        state0 = zeros(L, dtype=np.int)
        state1 = zeros(L, dtype=np.int)
        overflow0 = zeros(L, dtype=np.int)
        overflow1 = zeros(L, dtype=np.int)
        if len(init) == 2:
            state0[0] = init[0]
            state1[0] = init[1]

        # Implement the SDM
        for j in arange(1, L):
            state0[j] = state0[j - 1] + k[j - 1]
            if state0[j] > _maxvalue:
                overflow0[j] = 1
                state0[j] -= _maxvalue + 1
            state1[j] = state1[j - 1] + state0[j]
            if state1[j] > _maxvalue:
                overflow1[j] = 1
                state1[j] -= _maxvalue + 1
        sd = overflow0 + overflow1 - np.hstack(([0], overflow1[:-1]))
        state = state0 + state1
        cycles = np.where(state == 0)[0]
        if len(cycles) > 1:
            cycles = cycles[1] - cycles[0]
        else:
            cycles = -1

    # Modulator of order 3
    elif order == 3:
        # initialize the registers
        state0 = zeros(L, dtype=np.int)
        state1 = zeros(L, dtype=np.int)
        state2 = zeros(L, dtype=np.int)
        if len(init) == 3:
            state0[0] = init[0]
            state1[0] = init[1]
            state2[0] = init[2]
        overflow0 = zeros(L, dtype=np.int)
        overflow0[0] = 1
        overflow1 = zeros(L, dtype=np.int)
        overflow2 = zeros(L, dtype=np.int)

        # Implement the SDM
        for j in arange(1, L):
            state0[j] = state0[j - 1] + k[j - 1]
            if state0[j] > _maxvalue:
                overflow0[j] = 1
                state0[j] -= _maxvalue + 1

            state1[j] = state1[j - 1] + state0[j]
            if state1[j] > _maxvalue:
                overflow1[j] = 1
                state1[j] -= _maxvalue + 1

            state2[j] = state2[j - 1] + state1[j]
            if state2[j] > _maxvalue:
                overflow2[j] = 1
                state2[j] -= _maxvalue + 1
        sd = overflow0
        sd += overflow1 - np.hstack(([0], overflow1[:-1]))
        sd += overflow2 - 2 * np.hstack(([0], overflow2[:-1]))
        sd += np.hstack(([0], [0], overflow2[:-2]))

        state = state0 + state1 + state2
        cycles = np.where(state == 0)[0]
        if len(cycles) > 1:
            cycles = cycles[1] - cycles[0]
        else:
            cycles = -1

    elif order == 4:
        # initialize the registers
        state0 = zeros(L, dtype=np.int)
        state1 = zeros(L, dtype=np.int)
        state2 = zeros(L, dtype=np.int)
        state3 = zeros(L, dtype=np.int)
        if len(init) == 4:
            state0[0] = init[0]
            state1[0] = init[1]
            state2[0] = init[2]
            state3[0] = init[4]
        overflow0 = zeros(L, dtype=np.int)
        overflow0[0] = 1
        overflow1 = zeros(L, dtype=np.int)
        overflow2 = zeros(L, dtype=np.int)
        overflow3 = zeros(L, dtype=np.int)
        # Implement the SDM
        for j in arange(1, L):
            state0[j] = state0[j - 1] + k[j - 1]
            if state0[j] > _maxvalue:
                overflow0[j] = 1
                state0[j] -= _maxvalue + 1

            state1[j] = state1[j - 1] + state0[j]
            if state1[j] > _maxvalue:
                overflow1[j] = 1
                state1[j] -= _maxvalue + 1

            state2[j] = state2[j - 1] + state1[j]
            if state2[j] > _maxvalue:
                overflow2[j] = 1
                state2[j] -= _maxvalue + 1
            state3[j] = state3[j - 1] + state2[j]
            if state3[j] > _maxvalue:
                overflow3[j] = 1
                state3[j] -= _maxvalue + 1

        sd = (overflow0 +
              overflow1 - np.hstack(([0], overflow1[:-1])) +
              overflow2 - 2 * np.hstack(([0], overflow2[:-1])) +
              np.hstack(([0], [0], overflow2[:-2])) +
              overflow3 - 3 * np.hstack(([0], overflow3[:-1])) +
              3 * np.hstack(([0], [0], overflow3[:-2])) -
              np.hstack(([0], [0], [0], overflow3[:-3]))
              )

        state = state0 + state1 + state2 + state3
        cycles = np.where(state == 0)[0]
        if len(cycles) > 1:
            cycles = cycles[1] - cycles[0]
        else:
            cycles = -1

    return sd, cycles


def L_mash_dB(m, fref, n=1.0):
    """ Phase noise theoretical value of noise produced by a mash111 SDM

    This procedure calculates the noise at the output of the SD modulator

    Parameters
    ----------
    m : int
        The order of the SD modulator
    fref : float
        Reference frequency that is used to compare the output of the SD
        modulator
    n: float
        It is the average division ratio.

    return
    ------
    ldbc : Pnoise
        Return a function object of the noise
    """
    func_ldbc = lambda fm : 10 * log10((2 * pi) ** 2 / (12 * fref) * 
        (2 * sin(pi * fm / fref)) ** (2 * (m - 1)) / n ** 2)
    ldbc = Pnoise.with_function(func_ldbc, label='sdm')
    return ldbc 


if __name__ == "__main__":
    from numpy.testing import assert_almost_equal
    import numpy.random as rnd
    # Test order one assert the mean value
    floatnum = rnd.rand() * np.ones(100000)

    # order one
    sequence, cycles = gen_mash(1, 19, (floatnum * 2 ** 19).astype(int))
    assert_almost_equal(sequence.mean(), floatnum.mean(), 4)
    # order two
    sequence, cycles = gen_mash(2, 19, (floatnum * 2 ** 19).astype(int))
    assert_almost_equal(sequence.mean(), floatnum.mean(), 4)

    # order three
    sequence, cycles = gen_mash(3, 19, (floatnum * 2 ** 19).astype(int))
    assert_almost_equal(sequence.mean(), floatnum.mean(), 4)

    # order three
    sequence, cycles = gen_mash(3, 19, 0.25 * np.ones(100000) * 2 ** 19)
    assert_almost_equal(sequence.mean(), 0.25, 4)

    # order three
    sequence, cycles = gen_mash(4, 19, 0.25 * np.ones(100000) * 2 ** 19)
    assert_almost_equal(sequence.mean(), 0.25, 4)

    # Using the class the test is done as:
    sd_mash = SDModulator('mash', 3, 19,
                          (np.array([0.323232] * 100000) * 2 ** 19).astype(int))
    assert_almost_equal(sd_mash.seq.mean(), 0.323232, 4)
