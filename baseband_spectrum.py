import numpy as np
from numpy import  ceil, ones, linspace, cos, pi, arange, sqrt, shape, array
from numpy.fft import fft, fftshift
from scipy.signal import lfilter

def baseband_spectrum(x, Ts, BW, winval='blackman'):
    ''' function [Spectro, AF, B] = baseband_spectrum(x, Ts, BW):
        (Spectrum, AF, B) = my_baseband_spectrum(x, Ts, BW); estimates the spectrum 
     of X using the periodogram. The result is given in dBV (power integrated 
            % over BW). 
            
     It uses the Blackman window with an overlap of 50%.
    
     Windows
    
                   Main-lobe width         Side-lobe level
    
     Rect.         4*pi/N                  -13dB
    
     Barlett       8*pi/N                  -25dB
    
     Hamming       8*pi/N                  -41dB
    
     Blackman      12*pi/N                 -57dB''' 


    def rectangular():
        Aw_k = 2
        N = ceil(Aw_k/(BW*Ts))
        n = arange(0,N)
        w = ones(N)
        return(Aw_k,N,n,w)
    def barlett():
        Aw_k = 4
        N = ceil(Aw_k/(BW*Ts))
        n = arange(0,N)
        w = 1-abs((n+1-0.5*N)/N*0.5)
        return(Aw_k,N,n,w)
    def hamming():
        Aw_k = 4
        N = ceil(Aw_k/(BW*Ts))
        n = arange(0,N)
        w = 1-abs((t-0.5*N)/N*0.5)
        return(Aw_k,N,n,w)
    def blackman():
        Aw_k = 6;
        N = ceil(Aw_k/(BW*Ts))
        n = arange(1,N+1)
        w = 0.42 - 0.5*cos(2*pi*n/(N-1))+0.08*cos(4*pi*n/(N-1))
        return(Aw_k,N,n,w)
    window = {  'blackman' : blackman,
                'hamming'  : hamming,
                'rectangular': rectangular,
                'barlett': barlett,
            }
    # Fixed internal paramters
    winval = 'blackman'
    Aw_k,N,n,w =  window[winval]()
    overlap = 0.5;
    gas_pedal =  False;  
    # Scale the window by the mean squared value
    w = w/sqrt(np.dot(w,w)/len(w));
    
    L = shape(x);
    if N > L:
        raise Exception('baseband_spectrum:warning',
            'current number of samples'+ str(L)+'not enought.'+
            'Increase number of input samples to'+str(N))
    
    # Round the number of samples
    if (gas_pedal):
        N = int(2**(round(log2(N))))

    # Adjust periodgram method parameters
    AF = Aw_k/(N*Ts)  # Actual resolution bandwidth
    No = np.floor(N*overlap) # Number of overlapped samples
    B =  np.floor((len(x) - No)/(N - No));   # Number of dft's averaged
    # Spectrum calculation using the overlapped windowed periodogram
    aux = np.zeros(N)
    wx = w*x[0:N]
    aux = aux + abs(fft(wx))**2/N**2
    
    for b in np.arange(2,B+1):
        wx = w*x[(b-1)*(N-No):(b-1)*(N-No)+N];
        # wx  = wx - mean(wx); % Remove DC level in each frame
        aux = aux + abs(fft(wx))**2/N**2
    spectro = aux/B

    ## Normalization and averaging

    Aver = round(AF*N*Ts);   # Number of samples to perform averaging
    spectro = fftshift(spectro)
    # Averaging 
    spectro  = lfilter(ones(Aver), 1, spectro)
    # get rid of the edges
    spectro = np.r_[spectro[round(Aver/2)-1::],spectro[0:round(Aver/2)-1]]
    spectro[:round(Aver/2)] = spectro[round(Aver/2):2*round(Aver/2)]
    spectro[-round(Aver/2)-1:] = spectro[-2*round(Aver/2)-2:-round(Aver/2)-1]
    return(spectro,AF,B)

def main():
    Ts = 0.0001
    t = np.linspace(0,1,1/Ts);
    y = sin(2*np.pi*50*t) + sin(2*pi*120*t);
    y += x + 2*np.random.normal(scale=1, size=len(t));
    baseband_spectrum(y,Ts,1e3);
    
if __name__ == "__main__":
    main()
