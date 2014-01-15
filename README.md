plldesigner
===========

A pythonic tool for PLL design and exploration (focused in PLLs implemented in hardware). The final propose of this project is to have a complete design tool for PLL (Phase-locked loops) by creating a class that permits definer the different elements this class would eventually implement methods to:
* Analise the loop stability (that should be easy with scipy.signal)
* Specify the noise sources and calculate the overall noise
  - Using PWL file (trivial)
  - using the noise at 1Hz for the different components $1/f^0$,$1/f$,$1/f^2$,$1/f^3$ (not so usefull)
  - Integrate the noise using Gardner expression (done)
  - Import them from CVS or other type or format
  - Estimate the overall phase noise (by using a linear model)
  - Extend the plot utility to show diffent components and noise and the integrated noise
* Specify non-linearities in the loop and simulate the transient response (VCO, PFD-CP)
  - That is similar to what already is done by cppsim but this is a propertary software and I a want to see if is posible to create enough momentum to have such a simulator with a open source project.
  
Status
======

Currently this is in a embrionic state, for the first two bullets I understand what has to be done and I have a relativly complete software in Matlab that I have developed. Regarding the third bullet I have been thinking in different alternatives: generate veriloga code, generate modelica code, generate pure C++ or Cython code.


Develoment plant
================
The first answer to ask is which units to use. I think that is better to standarize everything externally to Single side banded noise (dBc/Hz),  this unit gives more insight that $\phi(f_m)$ or $\phi^2_(f_m)$

* Create the class structure to specify the noise,  the pnoise object, there have to be different init elements:
  1. Noise, specially for oscilator can be specified as $f_n/fm^{-n}+f_{n=1}/fm^{-n+1}+$ (done)
  2. over rule the addition operatior
