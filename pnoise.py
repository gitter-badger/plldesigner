# -*- coding: utf-8 -*-
"""
A Class to work with phase noise frequency data points
"""

from __future__ import (absolute_import, division, print_function)
import numpy as np
from numpy import log10, sqrt, sum
import matplotlib.pyplot as plt
import scipy.interpolate as intp


# functions to handle the units
__funits__ = {
    "dBc/Hz": lambda x: x,
    "rad/sqrt(Hz)": lambda x: 10 * log10(x ** 2 / 2),
    "rad**/Hz": lambda x: 10 * log10(x / 2),
}


class Pnoise(object):
    """
    Phase noise class to manipulate phase noise data in the frequency offset
    format

    Parameters
    ----------

    Atributes
    ---------


    """

    def __init__(self, fm, pnfm, label=None, units='dBc/Hz'):
        """
        Phase noise from separate vectors for fm and phase noise

        Parameters
        ----------
        fm : array_like
            Offset frequency vector
        pn : array_like
            Phase noise values

        Returns
        -------
        pnoise : Pnoise
        """
        self.units = units
        self.label = label
        self.fm = np.array(fm)

        # values for point slope approximation
        self.fi = None
        self.ldbc_fi = None
        self.slopes = None
        self.func_ldbc = None
        self.phi_out = None

        self.ldbc = __funits__[units](np.array(pnfm))

    @classmethod
    def with_interpolation(cls, fm, pnfm, label=None, units='dBc/Hz'):
        pnoise_class = cls(None, None, label=label, units='dBc/Hz')
        pnoise_class.fm = fm
        pnoise_class.ldbc = __funits__[units](np.array(pnfm))
        try:
            pnoise_class.func_ldbc = lambda fx: __pnoise_interp1d__(
                fm, pnoise_class.ldbc, fx)
        except ValueError:
            raise ValueError('The function can not be interpolated')
        return pnoise_class

    @classmethod
    def with_function(cls, func, label=None):
        """
        Class method for manipulating phase noise values

        Parameters
        ----------
        func : function
            A Python function representing the phase noise and function of
            the frequency offset: ldbc = func(fm) (dBc/Hz)
        label : string
            Name of the phase noise value

        Returns
        -------
        pnoise : Pnoise

        """
        pnoise_class = cls(None, None, label=label, units='dBc/Hz')
        pnoise_class.func_ldbc = func
        pnoise_class.label = label
        return pnoise_class

    @classmethod
    def with_points_slopes(cls, fi, ldbc_fi, slopes, label=None):
        """
        Phase noise with point and slope

        Parameters
        ----------
        fi : array_like
            Array with the offsets frequencies of the phase noise values
        ldbc_fi : array_like
            Array with the phase noise values at the fi frequencies
        slopes : array_like
            Array with slopes of the values that are interpolated (dBc/dec)
        label : str
            Plotting label

        Returns
        -------
        pnoise : Pnoise

        """
        pnoise_class = cls(fi, ldbc_fi, label=label, units='dBc/Hz')
        pnoise_class.slopes = slopes
        pnoise_class.func_ldbc = lambda fm: __pnoise_point_slopes__(
            fi, ldbc_fi, slopes, fm)
        return pnoise_class

    def eval_func(self, fx):
        pnoise_class = Pnoise(fx, self.func_ldbc(fx), label=self.label)
        return pnoise_class

    def plot(self, *args, **kwargs):
        plt.ylabel('$\mathcal{L}$(dBc/Hz)')
        plt.xlabel('$f_m$(Hz)')
        ax = plt.semilogx(self.fm, self.ldbc, label=self.label, *args,
                          **kwargs)
        return ax

    def __add__(self, other):
        """ Addition of though pnoise components """
        try:
            phi2fm = 2 * 10 ** (self.ldbc / 10)
            phi2fm_other = 2 * 10 ** (other.ldbc / 10)
            ldbc_add = 10 * log10((phi2fm + phi2fm_other) / 2)
        except ValueError:
            raise ValueError(
                'Additions is only allowed with vector of equal size'
                )
        add_noise = Pnoise(self.fm, ldbc_add)
        return add_noise

    def __mul__(self, mult):
        """ Multiplication of noise by a constant """
        if type(mult) not in (int, float, np.ndarray):
            raise TypeError('unsupported operand type(s) for mult')
        else:
            if type(mult) in (int, float):
                mult_noise = Pnoise(
                    self.fm, self.ldbc + 10 * log10(mult), label=self.label)
            else:
                try:
                    mult_noise = Pnoise(
                        self.fm, self.ldbc + 10 * log10(mult),
                        label=self.label)
                except ValueError:
                    raise ValueError('Vectors are not of the same length')
            return mult_noise

    def integrate(self, fl=None, fh=None, method='trapz'):
        """ Returns the integrated phase noise in (rad/rms)

            Parameters
            ----------
            fl :
            fh :
            method :

            Returns
            -------
            phi_out :

        """

        def gardner(ldbc_ix, fm_ix):
            """ Gardner integration method

                Parameters
                ----------
                ldbc_ix :
                fm_ix :

                Returns
                -------
                phi_sqrt :

            """
            lfm = len(ldbc_ix)
            # calculate the slope
            ai = ((ldbc_ix[1:lfm] - ldbc_ix[:lfm - 1]) /
                  (log10(fm_ix[1:lfm]) - log10(fm_ix[:lfm - 1])))
            """ If the slopes are never too big used Gardner method
            In simulations this is not the case """
            bi = (2 * 10 ** (ldbc_ix[:lfm - 1] / 10) * fm_ix[:lfm - 1] **
                  (-ai / 10) /
                  (ai / 10 + 1) * (fm_ix[1:lfm] ** (ai / 10 + 1) -
                                   fm_ix[:lfm - 1] ** (ai / 10 + 1)))
            if ai < 6:
                raise Warning(
                    'Gardner method fail due to large slope changes :\n' +
                    'use trapz method'
                )
            return sqrt(sum(bi))

        def trapz(ldbc_ix, fm_ix):
            """ Trapezoidal integration of the phase noise

                Parameters
                ----------
                ldbc_ix :
                fm_ix :

                Returns
                -------
                phi_out :

            """
            phi_2 = 2 * 10 ** (ldbc_ix / 10)
            return sqrt(np.trapz(phi_2, fm_ix))

        if fl is None:
            fl = min(self.fm)
        if fh is None:
            fh = max(self.fm)
        ix = (self.fm >= fl) & (self.fm <= fh)
        fi = self.fm[ix]
        ldbc_fi = self.ldbc[ix]
        if method == 'trapz':
            self.phi_out = trapz(ldbc_fi, fi)
        elif method == 'Gardner':
            self.phi_out = gardner(ldbc_fi, fi)
        else:
            raise ValueError('Integrating method not implemented')
        return self.phi_out


def __pnoise_interp1d__(fi, ldbc_fi, fm):
    """ Interpolate the phase noise assuming logaritmic linear behavior

        Parameters
        ---------
        fi :

        Returns
        ----------
        ldbc :

    """
    func_intp = intp.interp1d(log10(fi), ldbc_fi, kind='linear')
    ldbc = func_intp(log10(fm))
    return ldbc


def __pnoise_point_slopes__(fi, ldbc_fi, slopes, fm):
    """
    Function to evaluate a asymptotic model of the phase noise

    Parameters
    ----------
    fi :
    slopes :
    ldbc_fi :
    fm :

    returns
    -------
    ldbc_fm :

    """
    phi2 = 2 * 10 ** (ldbc_fi / 10)
    phi2 = phi2.reshape((1, len(phi2)))
    phi2_matrix = np.repeat(phi2, len(fm), axis=0)

    slopes = np.copy(slopes / 10)
    slopes = slopes.reshape((1, len(slopes)))
    slopes_matrix = np.repeat(slopes, len(fm), axis=0)

    fi = fi.reshape((1, len(fi)))
    fi_matrix = np.repeat(fi, len(fm), axis=0)

    fm = fm.reshape((len(fm), 1))
    fm_matrix = np.repeat(fm, fi.shape[1], axis=1)

    phi2_fm = np.sum(
        phi2_matrix * (fm_matrix / fi_matrix) ** slopes_matrix, axis=1)
    ldbc_fm = 10 * log10(phi2_fm / 2)
    return ldbc_fm


"""
Testing and other functions
"""
from numpy.testing import assert_almost_equal


def test__init__(plot=False):
    # Test __init__
    fm = np.logspace(3, 9, 100)
    lorentzian = Pnoise(fm, 10 * np.log10(1 / (fm * fm)), label='Lorentzian')
    white = Pnoise(fm, -120 * np.ones(fm.shape), label='white')
    added = white + lorentzian
    added.label = "addition"
    assert_almost_equal(added.ldbc[0], -60, 4)
    assert_almost_equal(added.ldbc[-1], -120, 4)
    ix, = np.where(fm > 1e6)
    assert_almost_equal(added.ldbc[ix[0]], -117.2822, 4)
    if plot:
        lorentzian.plot()
        white.plot()
        added.plot()
        plt.legend()
        plt.show()


def test_with_interpolation(plot=False):
    fm = np.logspace(3, 9, 1000)
    lorentzian = Pnoise.with_interpolation(
        fm, 10 * np.log10(1 / (fm * fm)), label='Lorentzian')
    fi = np.logspace(3, 9, 4)
    interpolated = lorentzian.eval_func(fi)
    ldbc_equation = 10 * np.log10(1 / (fi * fi))
    assert_almost_equal(ldbc_equation, interpolated.ldbc, 4)
    if plot:
        lorentzian.plot()
        interpolated.plot('o')
        plt.show()


def test_private_functions():
    # test the new
    fi = np.array([1e4, 1e9])
    ldbc_fi = np.array([-40, -150])
    slopes = np.array([-30, -20])
    fm = np.logspace(3, 9, 20)
    ldbc_model = __pnoise_point_slopes__(fi, ldbc_fi, slopes, fm)
    func = intp.interp1d(log10(fm), ldbc_model, kind='linear')
    ldbc_0 = func(log10(fi[0]))
    ldbc_1 = func(log10(fi[1]))
    assert_almost_equal(ldbc_0, ldbc_fi[0], 0)
    assert_almost_equal(ldbc_1, ldbc_fi[1], 0)


def test_with_points_slopes(plot=False):
    # test the new
    fi = np.array([1e4, 1e9])
    ldbc_fi = np.array([-40, -150])
    slopes = np.array([-30, -20])
    pnoise_model = Pnoise.with_points_slopes(fi, ldbc_fi, slopes)
    fm = np.logspace(3, 9, 20)
    pnoise_extrapolated = pnoise_model.eval_func(fm)
    if plot:
        pnoise_model.plot('o')
        pnoise_extrapolated.plot()
        plt.show()


def test_integration(plot=False):
    fm = np.logspace(4,8,1000)
    lorentzian  = Pnoise(fm,10*np.log10(1/(fm*fm)),label='Lorentzian')
    white = Pnoise(fm,-120*np.ones(fm.shape),label='white')
    added = lorentzian+white
    iadded_gardner = added.integrate()
    iadded_trapz = added.integrate(method='trapz')
    f_int = lambda fm : 2.0e-12*fm-2.0/fm
    i_f_int = np.sqrt(f_int(fm[-1])-f_int(fm[0]))
    assert_almost_equal(iadded_trapz, i_f_int, 5)
    assert_almost_equal(iadded_trapz, i_f_int, 5)


if __name__ == "__main__":
    test__init__(plot=False)
    test_private_functions()
    test_with_interpolation(plot=False)
    test_with_points_slopes(plot=False)
    test_integration(plot=False)