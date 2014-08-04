
import numpy as np
import matplotlib.pyplot as plt
import pll

class Sim(object):
    def __init__(self, resultdir, ts, tend):
        self.resuldir = resultdir
        self.ts = ts
        self.tend = tend 
        self.n = int(tend/ts)
        
class Struct(object):
    pass
        
if __name__ == "__main__":
    sim_tran = Sim('./sim/', 1e-9, 500e-6)
    phivco = np.zeros(sim_tran.n+1) + np.random.randn(1)
    phidiv = np.zeros(sim_tran.n+1)
    vcont =  np.zeros(sim_tran.n+1)
    fo =  np.zeros(sim_tran.n+1)
    vco = pll.vco(300e6, 26e9, 0.5)
    adpll = Struct()
    adpll.N = 57.23
    
    for i in np.arange(sim_tran.n)+1:
        fo[i] = vco.freq_vs_volt(vcont[i])
        phivco[i] = phivco[i-1]+fo[i]*sim_tran.ts
        phidiv[i] = phivco[i]/adpll.N
        
    # Remco tesis 
    # NsigmaNT = sim.Ts*VCO.Ffr;
    #sigmaNT  = sqrt(NsigmaNT * 10^(VCO.Lvco1M/10) * 1e6^2 / VCO.Ffr^3);
    # periodnoise = sigmaNT*randn(1,sim.n+1);
    plt.plot(phivco)
    plt.show()
    