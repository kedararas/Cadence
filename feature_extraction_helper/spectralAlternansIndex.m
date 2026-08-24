function [SAI, f_peak] = spectralAlternansIndex(metric)
%SPECTRALALTERNANSINDEX  Spectral alternans index of a per-beat metric series.
%
%   [SAI, f_peak] = spectralAlternansIndex(metric)
%
%   Power at the alternans frequency (0.4-0.5 cycles/beat) in the FFT of the
%   zero-mean, NaN-stripped per-beat series. Shared by analyzeAPAlternans
%   (Vm amplitude, APD80, dV/dt) and analyzeCaTransientAlternans (CaT
%   amplitude) so the Vm and Ca channels use one definition.
%
%   Returns SAI = 0, f_peak = 0.5 when fewer than 4 valid beats remain.

    m = metric(~isnan(metric));
    m = m - mean(m);
    N = length(m);
    if N < 4, SAI = 0; f_peak = 0.5; return; end

    Y    = fft(m);
    P    = abs(Y/N).^2;
    f    = (0:N-1)/N;
    band = f >= 0.40 & f <= 0.50;
    [SAI, i] = max(P(band));
    fb   = f(band);
    f_peak = fb(i);
end
