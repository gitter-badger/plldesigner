plldesigner
===========


A pythonic tool for PLL design and exploration (focused in PLLs implemented in hardware). More information can be found in [Phase Unlocked](http://jfosorio.github.io/). The final propose of this project is to have a complete design tool for PLL's (Phase-locked loops). It proposes a class that allows to:
* Analyze the loop stability 
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

Currently a first version is available The first two bullets are already cover. Regarding the third bullet I have been thinking in different alternatives: generate veriloga code, generate modelica code, generate pure C++ or Cython code.


Development plan
================


* Create the class structure to specify the noise,  the pnoise object, there have to be different init elements:
  1. Noise, specially for oscillator can be specified as $f_n/fm^{-n}+f_{n=1}/fm^{-n+1}+$ (done)
  2. over rule the addition operator (done)
  3. Interpolate (done)
  4. Generate a model of the data (extrapolate) (Pospone)
  5. Create plots with asymptotic values (Postpone)
  6. Plot several noise sources and the resultant 
* LTI model of the PLL
  1. Second order approximation (Done)
  2. Phase margin plot (To be Done)
  2. Timing vs phase  margin and error
  3. phase noise optimization 
* Use the design routines:
  1. Given fc and R (or the DN)  calculate the filter (done)
  2. Specify 
* Spectrum routines. (Those are not specifically needed only for this project but for others In any case can help in the implementation of the CP-PFD model) (I studied it
  pwell is enough) (Done)
* Circuit simulator
  1. First version with a fix step
     -CP is a big impulse as proposed by Perrot



