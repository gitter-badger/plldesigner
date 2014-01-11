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

port numpy as np
import scipy as sc
import matplotlib.pyplot as plt

class pnoise:
    ''' 
This is class defines objects and operations over phase noise values
There are different types of input functions:
    pwl
    ---
    pnoise(fm,fi,pn)
    pnoise(fm,fi,pn,function='pwl')
    pnoise(fm,pn,function='polynomic')
    with different units:
    pnoise(fm,fi,pn,
    ''' 
    def __init__(self,fm,*args,**kwargs):
        self.units = 'dBc/Hz'
        self.func = 'pwl'
        self.fm= np.array([])
        self.Ldbc=np.array([])
        if 'units' in kwargs:
            self.units = kwargs['units']
        if 'function' in kwargs:
            self.func = kawargs['function']
        self.func(fm,args,kwargs)

    def pwl(self,fm,pnfi):
        if self.units=='dBc/Hz':
            self.LdBc = pnfi 
            self.fm = fm
        elif self.units=='rad/sqrt(Hz)':
            self.fm = fm
            self.LdBc = 10*np.log10(pnfi**2/2)
        elif self.units=='rad**2/Hz':
            self.fm = fm
            self.LdBc = 10*np.log10(pnfi/2)

    def plot(self,*args,**kwargs):
        plt.ylabel('$\mathcal{L}$(dBc/Hz)')
        plt.xlabel('$f_m$(Hz)')
        fig= plt.semilogx(self.fm,10*np.log10(self.phi2fm/2),*args,**kwargs)
        return(fig)

    def __add__(self,other):
        add_noise = noise()
        add_noise.pwl(self.fm,self.phi2fm+other.phi2fm)
        return(add_noise)
'''
fm = np.logspace(3,9,100)
phi2 = noise(fm)
phi2.pwl(fm,fm*fm)
fig = phi2.plot('-o')
plt.show()
'''
