function results = mv_synth_cv(varargin)
%MV_SYNTH_CV  Validate the conduction-velocity engines against analytic truth.
%
%   results = mv_synth_cv
%   results = mv_synth_cv('Velocities', [20 60 100], 'Engines', {'polyfit'})
%
%   Conduction velocity is the one CADENCE metric that cannot be validated by
%   manual marking: a reviewer can mark activation times, but CV is a derived
%   spatial gradient with no hand-markable ground truth. This script supplies
%   that ground truth analytically, by building activation-time maps whose true
%   CV is known in closed form at every pixel, and comparing what
%   conduction_velocity() recovers.
%
%   The validation chain is therefore split in two, and both halves are covered:
%     * the activation-time INPUT  -> blinded manual marking (mv_mark)
%     * the gradient ALGORITHM     -> this script
%
%   TESTS
%     A  Planar wave, isotropic.  T = (x cos0 + y sin0)/v, so ||grad T|| = 1/v
%        exactly and true CV = v at every pixel, for any propagation angle.
%        Recovers speed and direction, and (polyfit) checks that the plane-fit
%        R^2 planarity metric reads ~1 on a genuinely planar wavefront.
%
%     B  Elliptical wave from a point source, anisotropic.  In fibre coordinates
%        T = sqrt((xi/vL)^2 + (eta/vT)^2), which gives longitudinal velocity vL
%        along the fibre axis, transverse vT across it, and a closed-form local
%        CV everywhere in between (see elliptic_T).  This is the test that
%        matters for the manuscript, because Figure 7 reports L, T and the
%        anisotropy ratio -- none of which a planar wave exercises.
%
%     C  Activation-time jitter sweep: Gaussian noise of known SD added to a
%        planar map, converting activation-marking precision (ms) into expected
%        CV error (%).  Gives the operating envelope and ties this script to the
%        manual-marking arm.
%
%     D  Cross-engine agreement: every test runs each engine on the SAME map, so
%        the summary doubles as a gradient-vs-polyfit comparison.
%
%   EXCLUSIONS (and why)
%     Edge pixels are excluded ('EdgeMargin'): the engines smooth the activation
%     map before differencing, which distorts the border.  Pixels near the point
%     source are excluded ('CoreRadius'): wavefront curvature there makes the
%     local CV genuinely lower, so counting them would score real physics as
%     algorithm error.  Both mirror what the CV module does on real data.
%
%   NOT COVERED HERE
%     The circle engine lives inside Cadence_Conduction_Velocity.mlapp
%     (circle_method_cv) rather than in conduction_velocity.m, so a script
%     cannot call it.  Extracting it to conduction_velocity_helper/ would let it
%     be tested alongside the other two.  Likewise the end-to-end test (synthetic
%     movie -> activation extraction -> CV) needs compute_lat_50, which is also
%     embedded in an .mlapp.
%
%   Name-value
%     'GridSize'    [nr nc] pixels (default [128 128])
%     'PixelSize'   mm per pixel (default 0.15)
%     'Velocities'  isotropic speeds to sweep, cm/s (default [20 40 60 80 100])
%     'Angles'      propagation angles, degrees (default [0 30 45 60 90])
%     'Aniso'       [vL vT] pairs, cm/s (default [60 25; 80 32; 40 20])
%     'Engines'     any of {'gradient','polyfit','circle'} (default all three)
%     'CircleRadius'  circle-engine radius in pixels (default 10)
%     'EdgeMargin'  border pixels excluded (default 12)
%     'CoreRadius'  pixels excluded around the point source (default 20)
%     'Smooth'      smooth_sigma_pix override; [] = engine default
%     'Tolerance'   pass threshold on median |error|, percent (default 5)
%     'Jitter'      activation-time jitter SDs to sweep, ms (default 0..4)
%     'JitterV'     planar speed used for the jitter sweep, cm/s (default 60)
%     'EndToEnd'    run the full-pipeline test E (default true)
%     'SNR'         SNRs for test E (default [Inf 20 10 5 3])
%     'SampleRate'  acquisition rate for test E, Hz (default 1000)
%     'PreSmooth'   spatial pre-smoothing sigmas for test E, px (default [0 2]),
%                   standing in for CADENCE's spatial denoising stage
%
%   Output
%     results  struct with .planar, .aniso and .jitter tables plus .passed

    p = inputParser;
    p.addParameter('GridSize',   [128 128], @(v) numel(v) == 2);
    p.addParameter('PixelSize',  0.15,      @isscalar);
    p.addParameter('Velocities', [20 40 60 80 100], @isnumeric);
    p.addParameter('Angles',     [0 30 45 60 90],   @isnumeric);
    p.addParameter('Aniso',      [60 25; 80 32; 40 20], @(v) size(v,2) == 2);
    p.addParameter('Engines',    {'gradient','polyfit','circle'}, @iscell);
    p.addParameter('CircleRadius', 10, @isscalar);
    p.addParameter('EdgeMargin', 12, @isscalar);
    p.addParameter('CoreRadius', 20, @isscalar);
    p.addParameter('Smooth',     [], @(v) isempty(v) || isscalar(v));
    p.addParameter('Tolerance',  5,  @isscalar);
    p.addParameter('Jitter',     [0 0.25 0.5 1 2 4], @isnumeric);
    p.addParameter('JitterV',    60, @isscalar);
    p.addParameter('EndToEnd',   true, @islogical);
    p.addParameter('SNR',        [Inf 20 10 5 3], @isnumeric);
    p.addParameter('SampleRate', 1000, @isscalar);
    p.addParameter('PreSmooth',  [0 2], @isnumeric);
    p.parse(varargin{:});
    o = p.Results;

    nr = o.GridSize(1);  nc = o.GridSize(2);  ds = o.PixelSize;

    fprintf('\n=========== mv_synth_cv ===========\n');
    fprintf('grid %dx%d px, pixel %.3f mm (field %.1f x %.1f mm)\n', ...
            nr, nc, ds, nc*ds, nr*ds);
    fprintf('engines: %s\n', strjoin(o.Engines, ', '));

    interior = false(nr, nc);
    interior(o.EdgeMargin+1:nr-o.EdgeMargin, o.EdgeMargin+1:nc-o.EdgeMargin) = true;

    % ================= TEST A: planar, isotropic =================
    fprintf('\n--- A: planar wave (isotropic), true CV = v everywhere ---\n');
    [eA, vA, aA, medA, p95A, dirA, r2A] = deal(strings(0,1), [], [], [], [], [], []);

    for e = 1:numel(o.Engines)
        eng = o.Engines{e};
        for iv = 1:numel(o.Velocities)
            v_cms = o.Velocities(iv);
            v_mm  = v_cms / 100;                     % cm/s -> mm/ms
            for ia = 1:numel(o.Angles)
                th = o.Angles(ia);
                T  = planar_T(nr, nc, ds, v_mm, th);

                [Vmag, Vx, Vy, ~, ~, meta] = run_engine(T, ds, eng, o.Smooth, o.CircleRadius);
                cv = Vmag * 100;                     % mm/ms -> cm/s

                m = interior & isfinite(cv);
                if ~any(m(:)), continue; end

                err  = 100 * (cv(m) - v_cms) / v_cms;
                dird = ang_err(atan2d(Vy(m), Vx(m)), th);

                eA(end+1,1)  = string(eng);          %#ok<AGROW>
                vA(end+1,1)  = v_cms;                %#ok<AGROW>
                aA(end+1,1)  = th;                   %#ok<AGROW>
                medA(end+1,1)= median(abs(err));     %#ok<AGROW>
                p95A(end+1,1)= prctile_(abs(err),95);%#ok<AGROW>
                dirA(end+1,1)= median(abs(dird));    %#ok<AGROW>
                if isfield(meta,'fit_r2') && ~isempty(meta.fit_r2)
                    r2A(end+1,1) = median(meta.fit_r2(m), 'omitnan'); %#ok<AGROW>
                else
                    r2A(end+1,1) = NaN;              %#ok<AGROW>
                end
            end
        end
    end

    planar = table(eA, vA, aA, medA, p95A, dirA, r2A, 'VariableNames', ...
        {'engine','v_true_cms','angle_deg','med_abs_err_pct','p95_abs_err_pct', ...
         'med_dir_err_deg','med_planarity_r2'});

    for e = 1:numel(o.Engines)
        s = planar.engine == string(o.Engines{e});
        if ~any(s), continue; end
        fprintf('  %-9s median |err| %.2f%%   p95 %.2f%%   direction %.2f deg', ...
                o.Engines{e}, median(planar.med_abs_err_pct(s)), ...
                max(planar.p95_abs_err_pct(s)), max(planar.med_dir_err_deg(s)));
        if any(isfinite(planar.med_planarity_r2(s)))
            fprintf('   planarity R2 %.4f', median(planar.med_planarity_r2(s),'omitnan'));
        end
        fprintf('\n');
    end

    % ================= TEST B: elliptical, anisotropic =================
    fprintf('\n--- B: elliptical wave (anisotropic), true CV known per pixel ---\n');
    [eB, vLB, vTB, medB, p95B, vLm, vTm, arT, arM] = ...
        deal(strings(0,1), [], [], [], [], [], [], [], []);

    r0 = round(nr/2);  c0 = round(nc/2);
    [cc, rr] = meshgrid(1:nc, 1:nr);
    core = hypot(cc - c0, rr - r0) <= o.CoreRadius;

    for e = 1:numel(o.Engines)
        eng = o.Engines{e};
        for k = 1:size(o.Aniso,1)
            vL_cms = o.Aniso(k,1);  vT_cms = o.Aniso(k,2);
            [T, cvTrue] = elliptic_T(nr, nc, ds, vL_cms/100, vT_cms/100, r0, c0);
            cvTrue = cvTrue * 100;                    % cm/s

            [Vmag, ~, ~, ~, ~, ~] = run_engine(T, ds, eng, o.Smooth, o.CircleRadius);
            cv = Vmag * 100;

            m = interior & ~core & isfinite(cv) & isfinite(cvTrue) & cvTrue > 0;
            if ~any(m(:)), continue; end
            err = 100 * (cv(m) - cvTrue(m)) ./ cvTrue(m);

            % Directional recovery along the two principal axes (fibre = +x).
            axL = false(nr,nc);  axL(r0, :) = true;   axL = axL & ~core & interior;
            axT = false(nr,nc);  axT(:, c0) = true;   axT = axT & ~core & interior;
            vLmeas = median(cv(axL & isfinite(cv)), 'omitnan');
            vTmeas = median(cv(axT & isfinite(cv)), 'omitnan');

            eB(end+1,1)  = string(eng);               %#ok<AGROW>
            vLB(end+1,1) = vL_cms;                    %#ok<AGROW>
            vTB(end+1,1) = vT_cms;                    %#ok<AGROW>
            medB(end+1,1)= median(abs(err));          %#ok<AGROW>
            p95B(end+1,1)= prctile_(abs(err),95);     %#ok<AGROW>
            vLm(end+1,1) = vLmeas;                    %#ok<AGROW>
            vTm(end+1,1) = vTmeas;                    %#ok<AGROW>
            arT(end+1,1) = vL_cms/vT_cms;             %#ok<AGROW>
            arM(end+1,1) = vLmeas/vTmeas;             %#ok<AGROW>
        end
    end

    aniso = table(eB, vLB, vTB, vLm, vTm, arT, arM, medB, p95B, 'VariableNames', ...
        {'engine','vL_true_cms','vT_true_cms','vL_meas_cms','vT_meas_cms', ...
         'ar_true','ar_meas','med_abs_err_pct','p95_abs_err_pct'});

    for i = 1:height(aniso)
        fprintf(['  %-9s  L %5.1f -> %5.1f   T %5.1f -> %5.1f   ' ...
                 'AR %.2f -> %.2f   field |err| med %.2f%%\n'], ...
                aniso.engine(i), aniso.vL_true_cms(i), aniso.vL_meas_cms(i), ...
                aniso.vT_true_cms(i), aniso.vT_meas_cms(i), ...
                aniso.ar_true(i), aniso.ar_meas(i), aniso.med_abs_err_pct(i));
    end

    % ================= TEST C: activation-time jitter =================
    % How far does error in the ACTIVATION INPUT propagate into CV? Gaussian
    % jitter is added to a planar map of known speed, so any deviation is purely
    % the gradient step amplifying LAT error. This is what links the two halves
    % of the validation: the manual arm measures how precisely activation can be
    % marked (an SD, in ms), and this table converts that SD into an expected CV
    % error. Noise is injected into the activation map rather than the optical
    % movie so that the CV algorithm is isolated from the upstream detector.
    fprintf('\n--- C: activation jitter -> CV error (planar, v = %g cm/s) ---\n', o.JitterV);
    [eC, jC, medC, p95C, keepC] = deal(strings(0,1), [], [], [], []);
    vj = o.JitterV / 100;
    Tclean = planar_T(nr, nc, ds, vj, 0);
    rng(11, 'twister');
    for e = 1:numel(o.Engines)
        eng = o.Engines{e};
        for k = 1:numel(o.Jitter)
            sd = o.Jitter(k);
            Tn = Tclean + sd * randn(nr, nc);
            [Vmag, ~, ~, ~, ~, ~] = run_engine(Tn, ds, eng, o.Smooth, o.CircleRadius);
            cv = Vmag * 100;
            m  = interior & isfinite(cv);
            if ~any(m(:))
                eC(end+1,1)=string(eng); jC(end+1,1)=sd; %#ok<AGROW>
                medC(end+1,1)=NaN; p95C(end+1,1)=NaN; keepC(end+1,1)=0; %#ok<AGROW>
                continue;
            end
            err = 100 * (cv(m) - o.JitterV) / o.JitterV;
            eC(end+1,1)   = string(eng);                        %#ok<AGROW>
            jC(end+1,1)   = sd;                                 %#ok<AGROW>
            medC(end+1,1) = median(abs(err));                   %#ok<AGROW>
            p95C(end+1,1) = prctile_(abs(err), 95);             %#ok<AGROW>
            keepC(end+1,1)= 100 * sum(m(:)) / sum(interior(:)); %#ok<AGROW>
        end
    end
    jitter = table(eC, jC, medC, p95C, keepC, 'VariableNames', ...
        {'engine','lat_jitter_ms','med_abs_err_pct','p95_abs_err_pct','pixels_kept_pct'});
    for i = 1:height(jitter)
        fprintf('  %-9s jitter %.2f ms -> med |err| %6.2f%%   p95 %7.2f%%   kept %5.1f%%\n', ...
                jitter.engine(i), jitter.lat_jitter_ms(i), jitter.med_abs_err_pct(i), ...
                jitter.p95_abs_err_pct(i), jitter.pixels_kept_pct(i));
    end

    % ================= TEST E: end-to-end pipeline =================
    % Everything above starts from an analytic activation map. This test starts
    % from a synthetic OPTICAL MOVIE, so the whole chain is exercised:
    %     AP waveforms -> compute_lat_50 -> conduction_velocity -> CV.
    %
    % A SINGLE beat is generated, because compute_lat_50 expects one windowed
    % beat (its search is confined to the depolarization phase of that window).
    % Feeding it a multi-beat trace lets max(dV/dt) latch onto the wrong beat.
    %
    % The sweep crosses SNR with a spatial PRE-SMOOTHING step standing in for
    % CADENCE's spatial denoising, because the interesting result is not the raw
    % noise tolerance but how sharply it depends on conditioning.
    %
    % UNITS: compute_lat_50 returns activation in FRAMES; conduction_velocity
    % documents T in ms. They coincide only at 1000 Hz, so LAT is converted
    % explicitly below (frames * 1000/fs).
    e2e = table();
    if o.EndToEnd
        fs = o.SampleRate;
        fprintf('\n--- E: end-to-end (movie -> compute_lat_50 -> CV), fs = %g Hz ---\n', fs);
        v_cms = o.JitterV;  v_mm = v_cms / 100;
        delay = planar_T(nr, nc, ds, v_mm, 0);
        fprintf('    adjacent-pixel LAT step %.3f ms; activation span %.1f ms\n', ...
                ds / v_mm, max(delay(:)));
        t_ms  = (0:round(200 * fs / 1000) - 1) / fs * 1000;
        base  = zeros(nr, nc, numel(t_ms));
        for rr2 = 1:nr
            for cc2 = 1:nc
                base(rr2,cc2,:) = ap_wave(t_ms, 20 + delay(rr2,cc2), 60);
            end
        end
        trueLat = 20 + delay;

        [eE, sE, pE, latsd, medE, p95E, keepE] = ...
            deal(strings(0,1), [], [], [], [], [], []);
        rng(23, 'twister');
        for k = 1:numel(o.SNR)
            snr = o.SNR(k);
            if isfinite(snr), stack0 = base + (1/snr)*randn(size(base));
            else,             stack0 = base; end
            for ps = o.PreSmooth
                stack = stack0;
                if ps > 0
                    for f = 1:size(stack,3)
                        stack(:,:,f) = imgaussfilt(stack(:,:,f), ps);
                    end
                end
                lat_ms = compute_lat_50(stack) * 1000 / fs;
                d = lat_ms(interior) - trueLat(interior);
                d = d(isfinite(d));
                sd = std(d - median(d));
                for e = 1:numel(o.Engines)
                    eng = o.Engines{e};
                    [Vmag, ~, ~, ~, ~, ~] = run_engine(lat_ms, ds, eng, o.Smooth, o.CircleRadius);
                    cv = Vmag * 100;
                    m  = interior & isfinite(cv);
                    eE(end+1,1) = string(eng); sE(end+1,1) = snr; pE(end+1,1) = ps; %#ok<AGROW>
                    latsd(end+1,1) = sd;                                            %#ok<AGROW>
                    if ~any(m(:))
                        medE(end+1,1)=NaN; p95E(end+1,1)=NaN; keepE(end+1,1)=0;     %#ok<AGROW>
                    else
                        err = 100 * (cv(m) - v_cms) / v_cms;
                        medE(end+1,1) = median(abs(err));                           %#ok<AGROW>
                        p95E(end+1,1) = prctile_(abs(err), 95);                     %#ok<AGROW>
                        keepE(end+1,1)= 100*sum(m(:))/sum(interior(:));             %#ok<AGROW>
                    end
                end
            end
        end
        e2e = table(eE, sE, pE, latsd, medE, p95E, keepE, 'VariableNames', ...
            {'engine','snr','presmooth_px','lat_sd_ms','med_abs_err_pct', ...
             'p95_abs_err_pct','pixels_kept_pct'});
        for i = 1:height(e2e)
            fprintf(['  %-9s SNR %5.1f  presmooth %gpx -> LAT SD %6.3f ms   ' ...
                     'CV |err| med %6.2f%%   kept %5.1f%%\n'], ...
                    e2e.engine(i), e2e.snr(i), e2e.presmooth_px(i), ...
                    e2e.lat_sd_ms(i), e2e.med_abs_err_pct(i), e2e.pixels_kept_pct(i));
        end
    end

    % ================= D: cross-engine agreement =================
    if numel(o.Engines) > 1
        fprintf('\n--- D: cross-engine (same maps, both engines) ---\n');
        for e = 1:numel(o.Engines)
            s = planar.engine == string(o.Engines{e});
            t = aniso.engine  == string(o.Engines{e});
            fprintf('  %-9s planar med |err| %.2f%%   anisotropic med |err| %.2f%%\n', ...
                    o.Engines{e}, median(planar.med_abs_err_pct(s)), ...
                    median(aniso.med_abs_err_pct(t)));
        end
    end

    % ================= verdict =================
    worst_planar = max(planar.med_abs_err_pct);
    worst_aniso  = max(aniso.med_abs_err_pct);
    passed = worst_planar <= o.Tolerance && worst_aniso <= o.Tolerance;

    fprintf('\n--- verdict (tolerance %.1f%% on median |error|) ---\n', o.Tolerance);
    fprintf('  worst planar      : %.2f%%\n', worst_planar);
    fprintf('  worst anisotropic : %.2f%%\n', worst_aniso);
    if passed
        fprintf('  PASS\n');
    else
        fprintf('  FAIL — an engine exceeds tolerance; inspect the tables\n');
    end
    fprintf('===================================\n');

    results = struct('planar', planar, 'aniso', aniso, 'jitter', jitter, ...
                     'e2e', e2e, 'passed', passed, 'opts', o);
end


% ======================================================================
function [Vmag, Vx, Vy, Gx, Gy, meta] = run_engine(T, ds, eng, sm, r)
%RUN_ENGINE  Dispatch to the requested CV engine with a common interface.
%   The circle engine lives in its own helper and takes a radius rather than
%   an options struct, so it is dispatched separately.
    if nargin < 5 || isempty(r), r = 10; end
    if strcmpi(eng, 'circle')
        [Vmag, Vx, Vy, Gx, Gy, meta] = circle_method_cv(T, ds, ds, r);
    else
        opt = struct('method', eng);
        if ~isempty(sm), opt.smooth_sigma_pix = sm; end
        [Vmag, Vx, Vy, Gx, Gy, meta] = conduction_velocity(T, ds, ds, opt);
    end
end


function T = planar_T(nr, nc, ds, v, theta)
%PLANAR_T  Activation map of a plane wave at speed v (mm/ms), angle theta (deg).
%   T = (x cos0 + y sin0)/v, so grad T has constant magnitude 1/v and the true
%   CV is exactly v at every pixel, independent of position and angle.
    [cc, rr] = meshgrid(0:nc-1, 0:nr-1);
    x = cc * ds;  y = rr * ds;                       % mm
    T = (x * cosd(theta) + y * sind(theta)) / v;     % ms
    T = T - min(T(:));                               % keep non-negative
end


function [T, cvTrue] = elliptic_T(nr, nc, ds, vL, vT, r0, c0)
%ELLIPTIC_T  Anisotropic activation map from a point source, fibres along +x.
%   T(xi,eta) = sqrt((xi/vL)^2 + (eta/vT)^2), the standard elliptical
%   approximation. Differentiating gives the closed-form local speed
%       CV = T / sqrt(xi^2/vL^4 + eta^2/vT^4)
%   which reduces to vL on the fibre axis (eta = 0) and vT across it (xi = 0),
%   and supplies per-pixel ground truth everywhere else. Undefined at the
%   source itself (0/0) — the caller excludes a core radius.
    [cc, rr] = meshgrid(1:nc, 1:nr);
    xi  = (cc - c0) * ds;                            % along fibre, mm
    eta = (rr - r0) * ds;                            % across fibre, mm
    T   = sqrt((xi / vL).^2 + (eta / vT).^2);        % ms
    den = sqrt(xi.^2 / vL^4 + eta.^2 / vT^4);
    cvTrue = T ./ max(den, eps);                     % mm/ms
end


function y = ap_wave(t_ms, ta, apd80)
%AP_WAVE  One action potential with APD80 exact by construction (see mv_selftest).
    t_up = ta + 2;
    tau  = (apd80 - 2) / log(5);
    y    = zeros(size(t_ms));
    rise = t_ms >= ta & t_ms < t_up;
    y(rise) = (t_ms(rise) - ta) / 2;
    dec = t_ms >= t_up;
    y(dec) = exp(-(t_ms(dec) - t_up) / tau);
end


function d = ang_err(a, b)
%ANG_ERR  Signed smallest angular difference, degrees, wrapped to [-180,180).
    d = mod(a - b + 180, 360) - 180;
end


function v = prctile_(x, p)
%PRCTILE_  Percentile without the Statistics Toolbox.
    x = sort(x(isfinite(x)));
    if isempty(x), v = NaN; return; end
    if isscalar(x), v = x; return; end
    pos = (p / 100) * numel(x) + 0.5;
    v   = interp1(1:numel(x), x, min(max(pos, 1), numel(x)), 'linear');
end
