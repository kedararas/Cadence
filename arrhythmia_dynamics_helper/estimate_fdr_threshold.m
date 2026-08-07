function [thr_frames, info] = estimate_fdr_threshold(lifespans, cycle_frames, target_fdr, fallback_frames)
%ESTIMATE_FDR_THRESHOLD  Data-driven PS-validity threshold at a target FDR.
%
%   [thr_frames, info] = estimate_fdr_threshold(lifespans, cycle_frames, ...
%                                               target_fdr, fallback_frames)
%
%   Chooses the smallest lifespan threshold (in frames) at which the estimated
%   false-discovery rate — the fraction of surviving tracks explainable by
%   chance-linked transient PS — drops to or below target_fdr.
%
%   Rationale: truly false PS (wavefront collisions, phase-map artifacts,
%   chance frame-to-frame linking) behave like a MEMORYLESS process, so their
%   track lifespans decay exponentially: S_noise(t) = exp(a - t/tau).  We fit
%   that decay to the short-lived (transient-dominated) regime and extrapolate.
%   At a candidate threshold the expected number of NOISE survivors is the
%   extrapolated S_noise; FDR = expected_noise / observed.  Real rotors decay
%   far slower, so observed >> expected_noise once t is past the noise tail.
%
%   Inputs
%     lifespans        vector of tracked-PS lifespans in frames (ps_info(:,3)).
%     cycle_frames     frames per rotation period (frame_rate / df_rotor).
%     target_fdr       desired max false-discovery rate (e.g. 0.10).
%     fallback_frames  threshold to use if the model can't be fit (e.g. the
%                      fixed min_cycles threshold).
%
%   Outputs
%     thr_frames   chosen threshold in frames.
%     info         struct: .tau, .a (fit params), .fdr_at_thr, .cycles_at_thr,
%                  .used_fallback (logical), .reason (text).
%
%   Notes
%     - A floor of 1 full rotation is enforced: a PS that has not completed one
%       rotation is not a rotor regardless of FDR arithmetic.
%     - With too few tracks or a degenerate fit the function returns
%       fallback_frames and flags info.used_fallback = true.

    lifespans = lifespans(:);
    lifespans = lifespans(isfinite(lifespans) & lifespans >= 0);

    info = struct('tau', NaN, 'a', NaN, 'fdr_at_thr', NaN, ...
                  'cycles_at_thr', NaN, 'used_fallback', false, 'reason', '');

    n_tracks = numel(lifespans);
    floor_frames = max(1, round(cycle_frames));   % >= 1 rotation

    % Need a reasonable sample to fit a decay curve.
    if n_tracks < 30 || ~isfinite(cycle_frames) || cycle_frames <= 0
        thr_frames = fallback_frames;
        info.used_fallback = true;
        info.reason = sprintf('too few tracks (%d) or invalid cycle length; used fallback', n_tracks);
        return;
    end

    % Survival curve S(t) = #tracks lasting >= t.
    t_max = max(lifespans);
    ts    = (1:t_max)';
    S     = arrayfun(@(t) sum(lifespans >= t), ts);

    % Fit the noise decay over the transient-dominated regime: 1 .. ~half a
    % rotation (short tracks are overwhelmingly false).  Require >= 6 points
    % with positive survival for a stable log-linear fit.
    t_fit_end = max(10, round(0.5 * cycle_frames));
    t_fit_end = min(t_fit_end, t_max);
    fit_idx   = ts <= t_fit_end & S > 0;
    if sum(fit_idx) < 6
        thr_frames = fallback_frames;
        info.used_fallback = true;
        info.reason = 'insufficient points in noise regime; used fallback';
        return;
    end

    % ln S = a - t/tau  ->  least squares for [a, 1/tau].
    A    = [ones(sum(fit_idx),1), -ts(fit_idx)];
    coef = A \ log(S(fit_idx));
    a      = coef(1);
    inv_tau = coef(2);
    if ~isfinite(inv_tau) || inv_tau <= 0
        thr_frames = fallback_frames;
        info.used_fallback = true;
        info.reason = 'non-decaying noise fit; used fallback';
        return;
    end
    tau = 1 / inv_tau;
    info.tau = tau;
    info.a   = a;

    % Walk thresholds from the 1-rotation floor upward; pick the first where
    % FDR <= target.
    chosen = NaN;
    for thr = floor_frames:t_max
        observed    = sum(lifespans >= thr);
        if observed == 0
            break;   % nothing survives beyond here
        end
        exp_noise = exp(a - thr / tau);
        fdr       = exp_noise / observed;
        if fdr <= target_fdr
            chosen = thr;
            info.fdr_at_thr    = fdr;
            info.cycles_at_thr = thr / cycle_frames;
            break;
        end
    end

    if isnan(chosen)
        % Target FDR never reached within the observed range — be conservative.
        thr_frames = max(fallback_frames, floor_frames);
        info.used_fallback = true;
        info.reason = 'target FDR not reached in observed range; used fallback';
        return;
    end

    thr_frames = chosen;
    info.reason = sprintf('FDR<=%.2f reached at %.2f cycles', target_fdr, info.cycles_at_thr);
end
