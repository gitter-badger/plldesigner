# -*- coding: utf-8 -*-
'''
First version of this class
Several open design questions:
    - in order to keep the same regular fm around
    - the class need to be interpolated and extrapolated
      how to do that ?
    - should I left the class just with the fm out of question ?
    - That would not allow to create noise from the corners
''' 

from __future__ import (absolute_import, division, print_function,
                                unicode_literals)

import numpy as np
from numpy import log10,sqrt,sum
import scipy as sc
import matplotlib.pyplot as plt

class pnoise:
    ''' 
This is class defines objects and operations over phase noise values
There are different types of input functions:
    pnoise(fm,LdBc)
    pnoise(fm,LdBc,label='label')
    pnoise(fm,phi,units='rad/sqrt(Hz)'
    The options of units are:
        dBc/Hz (default)
        rad/sqrt(Hz)
        rad**2/Hz
        ''' 
    def __init__(self,fm,pnfm,label=None,units='dBc/Hz'):
        self.units = units 
        self.label = label
        self.fm = fm 
        # functions to handle the units
        __funits = {
                'dBc/Hz' : lambda x: x,
                'rad/sqrt(Hz)' : lambda x: 10*log10(x**2/2), 
                'rad**/Hz' : lambda x: 10*log10(x/2), 
                }
        self.LdBc = __funits[units](pnfm)

    def plot(self,*args,**kwargs):
        plt.ylabel('$\mathcal{L}$(dBc/Hz)')
        plt.xlabel('$f_m$(Hz)')
        fig= plt.semilogx(self.fm,self.LdBc,label=self.label,*args,**kwargs)

    def __add__(self,other):
        try:
            phi2fm = 2*10**(self.LdBc/10)
            phi2fm_other = 2*10**(other.LdBc/10)
            LdBc_add = 10*log10((phi2fm+phi2fm_other)/2)
        except ValueError as er:
            print('Aditions is only allowed with vector of equal size')
        add_noise= pnoise(self.fm,LdBc_add)
        return(add_noise)

    def integrate(self,fl=[],fh=[]):
        '''Returns the integrated phase noise in rad over the limits fl,fh '''
        try:
            if fl==[]:
                fl = min(self.fm)
            if fh==[]:
                fh = max(self.fm)
            ix = (self.fm>=fl)  & (self.fm<=fh)
            fm_ix = self.fm[ix]
            LdBc_ix = self.LdBc[ix]
            lfm = len(LdBc_ix)
            ai =((LdBc_ix[1:lfm]-LdBc_ix[:lfm-1])/
                (log10(fm_ix[1:lfm])-log10(fm_ix[:lfm-1])))
            bi = (2*10**(LdBc_ix[:lfm-1]/10)*fm_ix[:lfm-1]**(-ai/10)/
                (ai/10+1)*(fm_ix[1:lfm]**(ai/10+1)-
                fm_ix[:lfm-1]**(ai/10+1)))
            self.phi_out = sqrt(sum(bi))
        except ValueError:
            print('verify fm and LBc have the same size')
        return(self.p)
        
'''
fm = nep.logspace(3,9,100)
lorentzian  = pnoise(fm,1/(fm*fm))
fig = lorentzian.plot('-o')
'''
