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
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

class noise:

    def __init__(self,fm,name="unamed",type='V^2/Hz'):
        ''' This is class defines objects and operations
        for noise objects'''
        self.type = type
        self.fm= np.array([])
        self.phi2fm=np.array([])

    def pwl(self,fm,phi):
        if self.type=='V/sqrt(Hz)':
            self.phi2fm = phi*phi 
            self.fm = fm
        else:
            self.fm = fm
            self.phi2fm = phi 

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
