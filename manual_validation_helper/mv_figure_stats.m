function stats = mv_figure_stats(metrics_file, varargin)
%MV_FIGURE_STATS  Numbers needed to lay out the arrhythmia figure.
%
%   stats = mv_figure_stats('/path/to/recording-metrics.mat')
%   stats = mv_figure_stats(metrics_file, 'Cam', 1, 'Channel', 'voltage')
%
%   Reports the three quantities that decide how the arrhythmia figure should
%   be built.  Each one answers a layout question that cannot be settled by
%   looking at the maps, because the thing being asked about is either below
%   the colour resolution of the image or hidden by the colour scale.
%
%   1. DOES THE RI MAP EARN A PANEL?
%      RI and OI share a denominator and differ only by harmonic content:
%      in cardiacSpectralMetrics/regOrgIndices,
%          RI = P[df +/- bw]                     / P[totBand]
%          OI = P[union of k*df +/- bw, k=1..K]  / P[totBand]
%      so OI >= RI always, and OI - RI IS the power at 2f, 3f, 4f.  A
%      near-sinusoidal rhythm (fast VF) has almost none, and the two maps are
%      then the same image twice.  This reports the OI - RI distribution over
%      tissue.  If the 95th percentile is small, RI carries no spatial
%      information the OI panel does not already show and belongs in the text
%      as a number, not in the figure as a map.
%
%   2. WHERE SHOULD THE DF COLOUR SCALE START AND STOP?
%      During entrained VF the whole field can sit in one FFT bin, which is a
%      real and reportable finding -- but on a wide colour scale it renders as
%      a flat silhouette and reads as an empty panel.  This reports the modal
%      DF, the fraction of tissue within one Welch bin of it (the homogeneity
%      claim, as a number worth quoting), and percentile-based limits that
%      show whatever spatial structure actually exists.
%
%   3. WHAT FRAME SPACING SHOWS EXACTLY ONE ROTATION?
%      A montage stepped in round numbers of milliseconds lands at an
%      arbitrary fraction of the cycle, so the last frame neither closes the
%      loop nor clearly fails to.  Spacing the frames at CL/(n-1) makes the
%      final frame reproduce the first, which is what demonstrates reentry
%      rather than merely sampling it.  This reports CL from the DF map and
%      the frame times and indices for a montage of n panels.
%
%   Reads only ep_metrics.complexity_data -- it never loads the camera stacks,
%   so it is fast on the multi-GB v7.3 metrics files.
%
%   Inputs
%     metrics_file  path to a *-metrics.mat written by Feature Extraction
%
%   Name-value
%     'Cam'       camera index (default 1)
%     'Channel'   'voltage' (default) or 'calcium'
%     'Frames'    panels in the rotor montage (default 5)
%     'AcqFreq'   sampling rate in Hz, only needed if the file does not
%                 carry acqFreq (default [] = read from file)
%     'RIThresh'  p95 of OI - RI below which RI is called redundant
%                 (default 0.10, i.e. RI within 10% of OI over 95% of tissue)
%
%   Output
%     stats  struct with fields .df, .ri, .oi, .montage, each holding the
%            summary values, plus .verdict_ri and .df_clim for the two
%            decisions.  Also printed to the console.
%
%   See also CARDIACSPECTRALMETRICS, MV_PACED_DF.

    p = inputParser;
    p.addParameter('Cam',      1,         @(v) isnumeric(v) && isscalar(v) && v >= 1);
    p.addParameter('Channel',  'voltage', @(s) any(strcmpi(s, {'voltage','calcium'})));
    p.addParameter('Frames',   5,         @(v) isnumeric(v) && isscalar(v) && v >= 2);
    p.addParameter('AcqFreq',  [],        @(v) isempty(v) || (isscalar(v) && v > 0));
    p.addParameter('RIThresh', 0.10,      @(v) isnumeric(v) && isscalar(v) && v > 0);
    p.parse(varargin{:});
    o = p.Results;

    %% ---------- load only the metrics struct ----------
    S = load(metrics_file, 'cmos_all_data');
    if ~isfield(S, 'cmos_all_data')
        error('mv_figure_stats:noData', ...
              '%s has no cmos_all_data struct.', metrics_file);
    end
    cad = S.cmos_all_data;

    if isempty(o.AcqFreq)
        if isfield(cad, 'acqFreq') && ~isempty(cad.acqFreq)
            fs = double(cad.acqFreq);
        else
            error('mv_figure_stats:noFs', ...
                  'acqFreq not in file; pass ''AcqFreq''.');
        end
    else
        fs = o.AcqFreq;
    end

    cx = cad.ep_metrics.complexity_data;
    if size(cx, 1) < o.Cam || isempty(cx{o.Cam, 1})
        error('mv_figure_stats:noCam', ...
              'No complexity_data for camera %d.', o.Cam);
    end

    % complexity_data{cam,1} is {DF, RI, OI}; on dual-channel recordings each
    % entry is itself {voltage, calcium}.
    maps = cx{o.Cam, 1};
    [DF, RI, OI] = unpack_maps(maps, o.Channel);

    tissue = isfinite(DF) & isfinite(RI) & isfinite(OI);
    n_tis  = nnz(tissue);
    if n_tis == 0
        error('mv_figure_stats:emptyMask', 'No finite pixels in the maps.');
    end

    %% ---------- 1. is RI redundant with OI? ----------
    % Same denominator by construction, so the difference is harmonic power.
    d   = OI(tissue) - RI(tissue);
    ri_stats = struct( ...
        'median_ri',   median(RI(tissue)), ...
        'median_oi',   median(OI(tissue)), ...
        'median_diff', median(d), ...
        'p95_diff',    prctile(d, 95), ...
        'max_diff',    max(d), ...
        'frac_neg',    mean(d < -1e-9));      % should be 0: OI >= RI always

    redundant   = ri_stats.p95_diff < o.RIThresh;
    verdict_ri  = ternary(redundant, ...
        'REDUNDANT - report RI as a number, do not give it a map panel', ...
        'INFORMATIVE - RI differs from OI, a map panel is justified');

    %% ---------- 2. DF homogeneity and colour limits ----------
    % Welch bin width matches cardiacSpectralMetrics' default NFFT.
    nfft     = max(256, 2^nextpow2(4*fs));
    bin_hz   = fs / nfft;

    dfv      = DF(tissue);
    modal_df = mode(round(dfv / bin_hz)) * bin_hz;   % snap to the bin grid
    frac_1bin = mean(abs(dfv - modal_df) <= bin_hz);

    % Percentile limits, widened to at least a few bins so a genuinely
    % single-bin map does not collapse to a zero-width colour axis.
    lo = prctile(dfv, 2);
    hi = prctile(dfv, 98);
    if hi - lo < 4*bin_hz
        mid = 0.5*(lo + hi);
        lo  = mid - 2*bin_hz;
        hi  = mid + 2*bin_hz;
    end

    df_stats = struct( ...
        'modal_df_hz',  modal_df, ...
        'median_df_hz', median(dfv), ...
        'bin_hz',       bin_hz, ...
        'frac_within_1bin', frac_1bin, ...
        'p2_hz',        prctile(dfv, 2), ...
        'p98_hz',       prctile(dfv, 98), ...
        'n_tissue_px',  n_tis);
    df_clim = [lo hi];

    %% ---------- 3. montage timing for one full rotation ----------
    % Use the modal DF: the montage shows the driver's rotation, and the
    % modal bin is what the field is entrained to.
    cl_ms       = 1000 / modal_df;
    step_ms     = cl_ms / (o.Frames - 1);       % last frame repeats the first
    times_ms    = (0:o.Frames-1) * step_ms;
    frame_idx   = round(times_ms * fs / 1000) + 1;

    montage = struct( ...
        'cl_ms',        cl_ms, ...
        'n_frames',     o.Frames, ...
        'step_ms',      step_ms, ...
        'times_ms',     times_ms, ...
        'frame_offset', frame_idx, ...
        'note', 'frame_offset are offsets from whichever frame you pick as t=0');

    %% ---------- pack and report ----------
    stats = struct('file', metrics_file, 'cam', o.Cam, 'channel', o.Channel, ...
                   'fs_hz', fs, 'df', df_stats, 'ri', ri_stats, ...
                   'oi_median', ri_stats.median_oi, 'df_clim', df_clim, ...
                   'montage', montage, 'verdict_ri', verdict_ri);

    print_report(stats, o);
end


% ========================================================================
function [DF, RI, OI] = unpack_maps(maps, channel)
%UNPACK_MAPS  Pull DF/RI/OI out of complexity_data, single- or dual-channel.
    if numel(maps) < 3
        error('mv_figure_stats:badCell', ...
              'complexity_data entry has %d elements, expected >= 3.', numel(maps));
    end
    pick = @(m) select_channel(m, channel);
    DF = pick(maps{1,1});
    RI = pick(maps{1,2});
    OI = pick(maps{1,3});
end


% ========================================================================
function m = select_channel(entry, channel)
%SELECT_CHANNEL  A dual-channel entry nests {voltage, calcium}; a
%   single-channel entry is already the map.
    if iscell(entry)
        idx = 1 + strcmpi(channel, 'calcium');
        if numel(entry) < idx
            error('mv_figure_stats:noChannel', ...
                  'No %s channel in this recording.', channel);
        end
        m = double(entry{idx});
    else
        m = double(entry);
    end
end


% ========================================================================
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end


% ========================================================================
function print_report(s, o)
    fprintf('\n=== mv_figure_stats: %s (CAM%d, %s) ===\n', ...
        s.file, s.cam, s.channel);
    fprintf('  sampling %.0f Hz, %d tissue pixels\n\n', s.fs_hz, s.df.n_tissue_px);

    fprintf('1) RI vs OI  (OI - RI = harmonic power; same denominator)\n');
    fprintf('     median RI          %.3f\n',  s.ri.median_ri);
    fprintf('     median OI          %.3f\n',  s.ri.median_oi);
    fprintf('     median OI-RI       %.3f\n',  s.ri.median_diff);
    fprintf('     p95    OI-RI       %.3f   (threshold %.2f)\n', ...
        s.ri.p95_diff, o.RIThresh);
    fprintf('     max    OI-RI       %.3f\n',  s.ri.max_diff);
    if s.ri.frac_neg > 0
        fprintf('     WARNING: %.2f%% of pixels have RI > OI, which the shared\n', ...
            100*s.ri.frac_neg);
        fprintf('              denominator makes impossible -- check the maps.\n');
    end
    fprintf('     VERDICT: %s\n\n', s.verdict_ri);

    fprintf('2) DF map\n');
    fprintf('     Welch bin          %.4f Hz\n', s.df.bin_hz);
    fprintf('     modal DF           %.3f Hz\n', s.df.modal_df_hz);
    fprintf('     median DF          %.3f Hz\n', s.df.median_df_hz);
    fprintf('     within 1 bin       %.1f%% of tissue\n', ...
        100*s.df.frac_within_1bin);
    fprintf('     2nd-98th pct       %.3f to %.3f Hz\n', s.df.p2_hz, s.df.p98_hz);
    fprintf('     suggested clim     [%.2f %.2f] Hz\n\n', ...
        s.df_clim(1), s.df_clim(2));

    fprintf('3) Rotor montage for one full rotation\n');
    fprintf('     cycle length       %.2f ms  (from modal DF)\n', s.montage.cl_ms);
    fprintf('     %d frames, step     %.2f ms\n', ...
        s.montage.n_frames, s.montage.step_ms);
    fprintf('     times (ms)         %s\n', ...
        strjoin(compose('%.1f', s.montage.times_ms), '  '));
    fprintf('     frame offsets      %s\n', ...
        strjoin(compose('%d', s.montage.frame_offset), '  '));
    fprintf('     (last frame reproduces the first -- that is the reentry check)\n\n');
end
