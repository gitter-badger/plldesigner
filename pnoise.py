# -*- coding: utf-8 -*-
'''
First version of this class
''' 

from __future__ import (absolute_import, division, print_function)

import numpy as np
from numpy import log10,sqrt,sum
import scipy as sc
import matplotlib.pyplot as plt
import scipy.interpolate as intp

class Pnoise(object):
    ''' 
This is class defines objects and operations over phase noise values
There are different types of input functions:
    Pnoise(fm,LdBc)
    Pnoise(fm,LdBc,label='label')
    Pnoise(fm,phi,units='rad/sqrt(Hz)'
    The options of units are:
        dBc/Hz (default)
        rad/sqrt(Hz)
        rad**2/Hz
        ''' 
    def __init__(self,fm,pnfm,label=None,units='dBc/Hz'):
        self.units = units 
        self.label = label
        self.fm = np.array(fm)
        # functions to handle the units
        __funits = {
                'dBc/Hz' : lambda x: x,
                'rad/sqrt(Hz)' : lambda x: 10*log10(x**2/2), 
                'rad**/Hz' : lambda x: 10*log10(x/2), 
                }
        self.LdBc = __funits[units](np.array(pnfm))

    def plot(self,*args,**kwargs):
        plt.ylabel('$\mathcal{L}$(dBc/Hz)')
        plt.xlabel('$f_m$(Hz)')
        ax= plt.semilogx(self.fm,self.LdBc,label=self.label,*args,**kwargs)
        return(ax)

    def __add__(self,other):
        try:
            phi2fm = 2*10**(self.LdBc/10)
            phi2fm_other = 2*10**(other.LdBc/10)
            LdBc_add = 10*log10((phi2fm+phi2fm_other)/2)
        except ValueError as er:
            print('Additions is only allowed with vector of equal size')
        add_noise= Pnoise(self.fm,LdBc_add)
        return(add_noise)
        
    def __mul__(self,mult):
        if type(mult) not in (int,float):
            raise TypeError('unsupported operand type(s) for mult')
        else:
            mult_noise = Pnoise(self.fm,self.LdBc+20*log10(mult))
            return(mult_noise)
            

    def interp1d(self,fi):
        '''Redifine the ordinate from the new fm to fi'''
        Lout = intp.interp1d(log10(self.fm),self.LdBc,log10(fi),'linear');
        pass

    def integrate(self,fl=[],fh=[],method='trapz'):
        '''Returns the integrated phase noise in rad over the limits fl,fh 
        Uses the Gardner algorithm.
        '''

        def gardner(LdBc_ix,fm_ix):
            ''' This is the Garder book integration method for the phase noise
            that does not work always with measurements or data'''
            lfm = len(LdBc_ix)
            #calculate the slope
            ai =((LdBc_ix[1:lfm]-LdBc_ix[:lfm-1])/
                (log10(fm_ix[1:lfm])-log10(fm_ix[:lfm-1])))
            if np.all(ai<6):
                """ If the slopes are never too big used Gardner method
                In simulations this is not the case """
                bi = (2*10**(LdBc_ix[:lfm-1]/10)*fm_ix[:lfm-1]**(-ai/10)/
                    (ai/10+1)*(fm_ix[1:lfm]**(ai/10+1)-
                    fm_ix[:lfm-1]**(ai/10+1)))
            return  sqrt(sum(bi))
             
        
        def trapz(LdBc_ix,fm_ix):
            phi_2 = 2*10**(LdBc_ix/10)
            return sqrt(np.trapz(phi_2,fm_ix))

        
        if fl==[]:
            fl = min(self.fm)
        if fh==[]:
            fh = max(self.fm)
        ix = (self.fm>=fl)  & (self.fm<=fh)
        fm_ix = self.fm[ix]
        LdBc_ix = self.LdBc[ix]
        if method=='trapz': 
            self.phi_out = trapz(LdBc_ix,fm_ix)
        elif method=='Gardner':
            self.phi_out = gardner(LdBc_ix,fm_ix)
        else: 
            raise Exception('Integrating method not implemented')
        return(self.phi_out)
        
'''
fm = nep.logspace(3,9,100)
lorentzian  = Pnoise(fm,1/(fm*fm))
fig = lorentzian.plot('-o')
'''
