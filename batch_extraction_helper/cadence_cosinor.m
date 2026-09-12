function C = cadence_cosinor(T, O, varargin)
%CADENCE_COSINOR  24 h cosinor test of every metric and onset measure across ZT.
%
%   C = cadence_cosinor(T, O)
%   C = cadence_cosinor(T, O, 'RefCL', 150, 'Condition', "Baseline", 'Period', 24)
%
%   For each metric in cadence_metric_labels the per-experiment value at the
%   reference CL (median over that experiment's recordings at that CL) is
%   fitted against ZT (h) with  y = M + A*cos(2*pi*t/P) + B*sin(2*pi*t/P).
%   The onset measures from cadence_alternans_onset (Onset_CL_V, Onset_CL_Ca,
%   Onset_CL_discordant, Arrhythmia_CL, FRP_CL, peak significant fractions)
%   are fitted the same way, one value per experiment.
%
%   Columns: Variable, RefCL, N, N_ZT, Mesor, Amplitude, Acrophase_h, R2, F, p,
%   Significant (p < 0.05).  Requires >= 6 experiments over >= 3 distinct ZTs.
%   p is the F-test of the cosinor model against a flat mean (2, N-3 d.f.).

    p = inputParser;
    p.addParameter('RefCL', 150);
    p.addParameter('Condition', "Baseline");
    p.addParameter('Period', 24);
    p.parse(varargin{:});
    o = p.Results;
    if ischar(T) || isstring(T), T = readtable(char(T), 'TextType', 'string', 'Delimiter', ','); end

    L = cadence_metric_labels();
    rows = {};

    cond = string(T.Condition); cond(ismissing(cond)) = "";
    sel = cond == string(o.Condition) & T.CL_ms == o.RefCL;
    key = strcat(string(T.ZT), "|", string(T.Experiment));
    for k = 1:size(L, 1)
        if ~ismember(L{k,1}, T.Properties.VariableNames), continue; end
        y = double(T.(L{k,1}));
        [t, v] = per_experiment(key(sel), T.ZT(sel), y(sel));
        rows(end+1, :) = fit_row(L{k,2}, o.RefCL, t, v, o.Period); %#ok<AGROW>
    end

    if nargin >= 2 && ~isempty(O) && istable(O)
        Osel = O(string(O.Condition) == string(o.Condition), :);
        onset_vars = {'Onset_CL_V','Onset_CL_Ca','Onset_CL_discordant','Onset_CL_concordant','Arrhythmia_CL', ...
                      'FRP_CL','Loss_of_capture_CL','Peak_V_sig_fraction','Peak_Ca_sig_fraction', ...
                      'Annotated_alternans_CL','Annotated_discordant_CL','Annotated_arrhythmia_CL'};
        for k = 1:numel(onset_vars)
            if ~ismember(onset_vars{k}, Osel.Properties.VariableNames), continue; end
            rows(end+1, :) = fit_row(['[onset] ' onset_vars{k}], NaN, double(Osel.ZT), double(Osel.(onset_vars{k})), o.Period); %#ok<AGROW>
        end
    end

    C = cell2table(rows, 'VariableNames', {'Variable','RefCL','N','N_ZT','Mesor','Amplitude','Acrophase_h','R2','F','p','Significant'});
end


function [t, v] = per_experiment(key, zt, y)
    [uk, ~, g] = unique(key, 'stable');
    t = nan(numel(uk), 1); v = t;
    for i = 1:numel(uk)
        idx = g == i;
        t(i) = median(zt(idx), 'omitnan');
        yy = y(idx); yy = yy(isfinite(yy));
        if ~isempty(yy), v(i) = median(yy); end
    end
end

function r = fit_row(name, refcl, t, y, P)
    ok = isfinite(t) & isfinite(y);
    t = t(ok); y = y(ok);
    n = numel(y); nzt = numel(unique(t));
    r = {name, refcl, n, nzt, NaN, NaN, NaN, NaN, NaN, NaN, ""};
    if n < 6 || nzt < 3, return; end
    if var(y) <= 1e-20 * max(1, mean(y)^2), return; end     % flat: nothing to fit
    X = [ones(n,1), cos(2*pi*t/P), sin(2*pi*t/P)];
    b = X \ y;
    yhat = X * b;
    sse = sum((y - yhat).^2); sst = sum((y - mean(y)).^2);
    r2 = 1 - sse / sst;
    df2 = n - 3;
    if df2 <= 0
        F = NaN; pval = NaN;
    elseif sse <= 1e-12 * sst
        F = Inf; pval = 0;                       % perfect fit
    else
        F = max(0, ((sst - sse) / 2) / (sse / df2));   % clamp rounding (sse ~ sst -> F ~ -eps)
        pval = 1 - fcdf_local(F, 2, df2);
    end
    amp = hypot(b(2), b(3));
    acro = mod(atan2(b(3), b(2)) / (2*pi) * P, P);
    r = {name, refcl, n, nzt, b(1), amp, acro, r2, F, pval, string(ifelse(pval < 0.05, "yes", "no"))};
end

function s = ifelse(c, a, b)
    if c, s = a; else, s = b; end
end

function p = fcdf_local(F, d1, d2)
% F cumulative distribution via the regularized incomplete beta (no toolbox needed)
    if ~isfinite(F), p = NaN; return; end
    x = min(max(d1 * F / (d1 * F + d2), 0), 1);
    p = betainc(x, d1/2, d2/2);
end
