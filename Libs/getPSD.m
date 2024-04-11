% Gets the Power Spectral Density (PSD) of the signal acc.
%
% Parameters:
%   - acc: the input signal.
%   - dt: the time step between samples.
% Returns:
%   - f: the frequency vector.
%   - PSD: the PSD of the signal acc.
function [f, PSD] = getPSD(acc, dt)

    % Remove the singleton dimensions
    acc = squeeze(acc);
    Fs = 1 / dt;
    L = length(acc);
    NFFT = 2 ^ nextpow2(L);

    % Calculate the PSD
    Y = fft(acc, NFFT) / L;
    f = Fs / 2 * linspace(0, 1, NFFT / 2 + 1);
    PSD = 2 * abs(Y(1:NFFT / 2 + 1));
