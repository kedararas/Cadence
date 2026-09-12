function [power_map, k_map, ok_map, nUsed] = spectralAlternansMap(beats, minBeatFrac)
%SPECTRALALTERNANSMAP  Per-pixel spectral alternans index of a per-beat stack.
%
%   [power_map, k_map, ok_map, nUsed] = spectralAlternansMap(beats)
%   [...] = spectralAlternansMap(beats, minBeatFrac)
%
%   Per-pixel counterpart of SPECTRALALTERNANSINDEX, shared by
%   analyzeAPAlternans and analyzeCaTransientAlternans so the Vm and Ca
%   channels use one definition.
%
%   INPUT
%     beats        R x C x nBeats stack of a per-beat metric (APD80, CaT
%                  amplitude, ...), NaN wherever that beat was not measurable
%                  at that pixel.
%     minBeatFrac  fraction of the used beats that must be finite for a pixel
%                  to contribute (default 1.0 = every beat).  1.0 is the safe
%                  setting: no substituted value enters the transform at all.
%                  Relaxing it keeps more tissue but readmits the substitution
%                  artefact described below, so lower it only while watching
%                  the coverage of ok_map.
%
%   OUTPUTS
%     power_map  power at 0.5 cycles/beat, NaN outside ok_map.  Same definition
%                and units as the value this function replaces, but masked.
%     k_map      dimensionless Nearing-Verrier style score: the alternans bin
%                minus the mean of the surrounding noise bins, over their SD.
%                Comparable across recordings and between Vm and Ca, which raw
%                power is not -- power scales with signal amplitude, which is
%                why non-zero sai_v ran ~2e-2 while ca_sai ran ~2e-6, a units
%                artefact rather than biology.
%     ok_map     logical mask of the pixels that contributed.
%     nUsed      beats actually transformed (even; see below).
%
%   WHY THIS EXISTS
%     The per-pixel map used to be built in-line as
%         spectral_map = nan(R, C);
%         all_x(isnan(all_x)) = 0;              % background NaN -> finite 0
%         ...
%         spectral_map = P(:,:, alt_bin);       % overwrites the whole map
%     which turned every out-of-tissue NaN into a finite zero and then replaced
%     the NaN-initialised map wholesale.  The pipeline contract is that maps
%     arrive masked and cadence_recording_medians' nanmed() simply drops
%     non-finite pixels -- so the reported median was taken over the entire
%     frame with roughly half of it hard zero.  Once background plus invalid
%     tissue passed 50% of the frame, the median landed exactly on 0.  Calcium
%     masks sit below half the frame (median tissue fraction 0.488, vs 0.560
%     for voltage), so ca_sai was exactly zero in ~79% of baseline recordings
%     -- 100% below 60 ms -- and sai_v in ~37%.  Every other map in these
%     functions escaped this because it is built from omitnan arithmetic on
%     NaN inputs, which propagates NaN.
%
%     Three further defects, also fixed here:
%       * NaN->0 on a pixel valid for only some beats injects a large
%         artificial excursion, and when beat validity itself alternates (2:1
%         capture, or alternate beats failing detection) that substitution
%         deposits power at exactly the alternans bin -- manufacturing the
%         signal being measured.  The mean is now taken over valid beats only,
%         and partial pixels are excluded by default.
%       * floor(nBeats/2)+1 is the exact 0.5 cycles/beat bin only when nBeats
%         is even; at 47 beats it selects 23/47 = 0.489 cycles/beat, so true
%         alternans power is split across neighbouring bins and under-read.
%         The stack is truncated to an even beat count, which is what alternans
%         -- defined on beat pairs -- wants in any case.
%       * raw power is reported alongside, but k_map is the comparable form.
%
%   See also SPECTRALALTERNANSINDEX, ANALYZEAPALTERNANS,
%   ANALYZECATRANSIENTALTERNANS.

    if nargin < 2 || isempty(minBeatFrac), minBeatFrac = 1.0; end

    [R, C, nBeats] = size(beats);
    power_map = nan(R, C);
    k_map     = nan(R, C);
    ok_map    = false(R, C);

    nUsed = nBeats - mod(nBeats, 2);          % even -> exact 0.5 cyc/beat bin
    if nUsed < 4, return; end

    A     = beats(:, :, 1:nUsed);
    fin   = isfinite(A);
    valid = sum(fin, 3);

    if minBeatFrac >= 1
        need = nUsed;
    else
        need = max(4, ceil(minBeatFrac * nUsed));
    end
    ok_map = valid >= need;
    if ~any(ok_map(:)), return; end

    A(~fin) = 0;
    % Mean over the VALID beats only.  A plain mean(A,3) after NaN->0 drags the
    % baseline of every partial pixel toward zero, and that is precisely what
    % injected false power at the alternans bin.
    A = A - (sum(A, 3) ./ max(valid, 1));
    A(~repmat(ok_map, 1, 1, nUsed)) = 0;

    P       = abs(fft(A, [], 3) / nUsed).^2;
    alt_bin = nUsed/2 + 1;                    % 1-based: bin k holds (k-1)/nUsed
    alt     = P(:, :, alt_bin);

    noise_bins = 2:(alt_bin - 1);             % exclude DC and the alternans bin
    if numel(noise_bins) >= 3
        nmu = mean(P(:, :, noise_bins), 3);
        nsd = std(P(:, :, noise_bins), 0, 3);
        k   = (alt - nmu) ./ (nsd + eps);
    else
        k = nan(R, C);                        % too few beats for a noise floor
    end

    alt(~ok_map) = NaN;
    k(~ok_map)   = NaN;
    power_map = alt;
    k_map     = k;
end
