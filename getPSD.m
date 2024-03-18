function [f,PSD] = getPSD(acc,dt)
    %Gets the PSD of the signal acc
    acc= squeeze(acc);
    Fs = 1/dt;
    L = length(acc);
    NFFT = 2^nextpow2(L);
    Y = fft(acc,NFFT)/L;
    f = Fs/2*linspace(0,1,NFFT/2+1);
    PSD = 2*abs(Y(1:NFFT/2+1));
