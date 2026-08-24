function [Vmag, Vx, Vy, Gx, Gy, meta] = conduction_velocity(T, dx, dy, opts)
% CONDUCTION_VELOCITY  Compute conduction velocity from an activation-time map.
%
% Inputs
%   T    : activation time map (2D, ms). NaN/0 = off-tissue (same convention
%          as circle_method_cv; masked maps use a 0 background).
%          Any consistent time unit works; the CADENCE app passes ms.
%   dx   : spatial step in x-direction (e.g., mm per pixel)
%   dy   : spatial step in y-direction (e.g., mm per pixel)
%   opts : struct with optional fields:
%          .method            : 'gradient' (finite-difference of smoothed T,
%                               default) or 'polyfit' (grad from a local least-
%                               squares surface fit — more robust, curvature-aware)
%          .smooth_sigma_pix  : Gaussian sigma in pixels (default 0 = none)
%          .min_support       : reject smoothed pixels whose valid fraction of
%                               kernel mass is below this (default 0.5; 1 =
%                               fully valid neighbourhood, ~0.5 = straight edge)
%          .grad_min          : minimum |grad T| to trust (ms/cm) (default auto)
%          .speed_max         : cap CV magnitude (mm/ms) (default 2.50)
%          .polyfit_win       : half-window (pixels) for 'polyfit' (default 3 -> 7x7)
%          .polyfit_order     : 1 = plane, 2 = quadratic (default 2, captures curvature)
%          .use_imgaussfilt   : true to use imgaussfilt if available (default true)
%
% Outputs
%   Vmag : CV magnitude (mm/ms == m/s; the caller multiplies by 100 for cm/s)
%   Vx   : CV x-component (mm/ms)
%   Vy   : CV y-component (mm/ms)
%   Gx   : dT/dx (ms/mm)
%   Gy   : dT/dy (ms/mm)
%   meta : struct with masks & settings (valid mask, thresholds, etc.)
%
% Relationship used:
%   speed  = 1 / ||∇T||,
%   vector = -∇T / ||∇T||^2
%   (direction points from early->late; units consistent when T is in ms and
%   x,y in mm, giving CV in mm/ms)

if nargin < 4, opts = struct(); end
opts = setdefault(opts, 'method', 'gradient');     % 'gradient' | 'polyfit'
opts = setdefault(opts, 'smooth_sigma_pix', 2);
opts = setdefault(opts, 'min_support', 0.5);       % min valid fraction of kernel mass (1 = fully valid, ~0.5 = straight edge)
opts = setdefault(opts, 'grad_min', []);
opts = setdefault(opts, 'speed_max', 2.50);        % adjust to your prep/species
opts = setdefault(opts, 'reject_over_cap', true);  % true: drop over-cap pixels; false: clamp to speed_max
opts = setdefault(opts, 'polyfit_win', 3);         % half-window (pixels) for 'polyfit'
opts = setdefault(opts, 'polyfit_order', 2);       % 1 = plane, 2 = quadratic
opts = setdefault(opts, 'use_imgaussfilt', true);

T = double(T);
% Shared tissue convention with circle_method_cv: NaN/0 = off-tissue. Masked
% activation maps are produced by multiplication (act_times .* mask), so a
% 0 background means "outside the mask", not "activated at t=0"; without this
% the gradient engine smooths real LATs into the background and emits
% spurious border vectors that the circle engine (correctly) rejects.
T(T == 0) = NaN;
validT = isfinite(T);

% --- Optional smoothing (on T, not on gradients) ---
% NaN-AWARE (normalised / Nadaraya-Watson convolution).  A plain Gaussian is
% NOT NaN-aware: every pixel whose kernel touches a NaN becomes NaN, so at the
% default sigma=2 the valid region erodes ~3*sigma = 6 px inward from every hole
% AND from the tissue boundary.  On a masked activation map that is the dominant
% cause of missing CV pixels -- far larger than the grad_min or speed_max
% rejections.  Measured on real data before this fix: a human wedge LAT map lost
% 38.6% of its finite pixels to the filter alone, and CV coverage fell from
% 97.5% (sigma=0) to 39.7% (sigma=2) with the loss being one large contiguous
% region, not speckle.
%
% Instead, filter the NaN-zeroed map and the validity mask with the SAME kernel
% and divide.  Each output pixel is then the Gaussian-weighted mean of the VALID
% neighbours only, so smoothing no longer eats coverage.  Near a boundary the
% kernel is one-sided, which is the standard (and far preferable) behaviour --
% callers that need edge pixels excluded should still apply their own margin.
if opts.smooth_sigma_pix > 0
    sig = opts.smooth_sigma_pix;
    Tf  = T;  Tf(~validT) = 0;      % zero-fill: contributes nothing to the sum
    W   = double(validT);           % ...and its weight is zero in the denominator

    if opts.use_imgaussfilt && exist('imgaussfilt','file')
        % Safe now: Tf and W contain no NaN, so imgaussfilt cannot spread one.
        num = imgaussfilt(Tf, sig, 'FilterDomain', 'spatial');
        den = imgaussfilt(W,  sig, 'FilterDomain', 'spatial');
    else
        % fallback: separable Gaussian via conv2
        ksz = max(3, 2*ceil(3*sig)+1);
        g = exp(-((-(ksz-1)/2:(ksz-1)/2).^2)/(2*sig^2));
        g = g / sum(g);
        num = conv2(conv2(Tf, g, 'same'), g', 'same');
        den = conv2(conv2(W,  g, 'same'), g', 'same');
    end

    Ts = num ./ den;
    % Support floor. den is the kernel-weighted VALID fraction of each
    % neighbourhood (kernels are normalised: 1 = fully valid, ~0.5 = straight
    % tissue edge). A bare den<=eps guard can never fire for a valid pixel —
    % its own kernel weight (~0.04 at sigma=2) already exceeds eps — so
    % near-isolated pixels kept essentially unsmoothed values and pixels whose
    % kernel reaches across NaN holes (e.g. block lines) blended times from
    % the far side. Pixels with less than min_support valid mass are rejected;
    % den is exported as meta.smooth_support so callers can gate stricter.
    Ts(den < opts.min_support) = NaN;
    support = den;
else
    Ts = T;
    support = [];
end

% Preserve NaN regions; support-rejected pixels are no longer valid either
% (keeps them out of polyfit windows and downstream masks).
Ts(~validT) = NaN;
validT = validT & isfinite(Ts);

% --- Gradients ---
% Both engines feed the same speed = 1/||grad T|| relationship below; they only
% differ in HOW grad T is estimated.
fitR2 = [];
if strcmpi(opts.method, 'polyfit')
    % Grad from a local least-squares polynomial surface fit over a window of
    % VALID neighbours (Bayly et al. 1998). Regularizes without the per-pixel
    % noise amplification of finite differences; order 2 also captures wavefront
    % curvature. fitR2 is the per-pixel goodness-of-fit (a quality metric).
    [Gx, Gy, fitR2] = local_polyfit_gradient(Ts, validT, dx, dy, ...
        opts.polyfit_win, opts.polyfit_order);
else
    % MATLAB's [FX,FY]=gradient(F,hx,hy) pairs the 1st output (x/columns) with hx
    % and the 2nd (y/rows) with hy, so pass dx then dy. (Currently dx==dy from the
    % app, but keeping the pairing correct guards against anisotropic pixels.)
    [Gx, Gy] = gradient(Ts, dx, dy);  % Gx=dT/dx (ms/mm), Gy=dT/dy (ms/mm)
end

% --- Gradient magnitude and validity ---
gradMag = hypot(Gx, Gy);          % ||∇T|| (ms/mm)
% Auto threshold if not provided: small fraction of robust central tendency
if isempty(opts.grad_min)
    gm = gradMag(validT);
    gm = gm(isfinite(gm));
    if isempty(gm)
        opts.grad_min = 0;
    else
        % 2.5% of median: the applied threshold. (Previously written as 5% of
        % median but applied at half strength; folded so the auto default is
        % unchanged while opts.grad_min is now applied — and reported in
        % meta — exactly as given.)
        opts.grad_min = max(5e-7, 0.025*median(gm));
    end
end

good = validT & isfinite(gradMag) & (gradMag >= opts.grad_min);

% --- Conduction velocity (mm/ms) ---
Vmag = nan(size(T));
Vx   = nan(size(T));
Vy   = nan(size(T));

Vmag(good) = 1 ./ gradMag(good);
% Vector form: v = -∇T / ||∇T||^2
%Vx(good)   = -Gx(good) ./ (gradMag(good).^2);
Vx(good)   = Gx(good) ./ (gradMag(good).^2);
Vy(good)   = Gy(good) ./ (gradMag(good).^2);
%Vy(good)   = -Gy(good) ./ (gradMag(good).^2);

% Handle unphysiological speeds. Over-cap pixels come from near-flat LAT
% (||grad T|| -> 0 gives 1/||grad T|| -> Inf), i.e. the wavefront plateau, not
% real fast conduction. Clamping them to speed_max leaves a plateau of fake
% max-speed vectors that pollute directional/anisotropy pooling and the quiver;
% rejecting them (default) marks those pixels invalid instead.
if ~isempty(opts.speed_max) && isfinite(opts.speed_max)
    tooFast = Vmag > opts.speed_max;
    if opts.reject_over_cap
        Vmag(tooFast) = NaN;
        Vx(tooFast)   = NaN;
        Vy(tooFast)   = NaN;
        good(tooFast) = false;
    else
        Vmag(tooFast) = opts.speed_max;
        scale = opts.speed_max ./ hypot(Vx, Vy);
        scale(~isfinite(scale)) = 1;
        Vx(tooFast) = Vx(tooFast) .* scale(tooFast);
        Vy(tooFast) = Vy(tooFast) .* scale(tooFast);
    end
end

meta = struct();
meta.method     = lower(opts.method);   % lets ensure_local_cv_field cache-key the engine
meta.validT     = validT;
meta.good       = good;
meta.grad_min   = opts.grad_min;
meta.speed_max  = opts.speed_max;
meta.reject_over_cap = opts.reject_over_cap;
meta.dx         = dx;
meta.dy         = dy;
meta.smooth_sigma_pix = opts.smooth_sigma_pix;
meta.min_support      = opts.min_support;
meta.smooth_support   = support;             % valid kernel-mass fraction per pixel ([] if no smoothing)
meta.fit_r2     = fitR2;                 % [] for gradient; per-pixel fit R^2 for polyfit

end

function [Gx, Gy, Rsq] = local_polyfit_gradient(T, valid, dx, dy, win, order)
% Per-pixel gradient from a local least-squares polynomial surface fit.
% Fits T(x,y) ~ poly(order) over the (2*win+1)^2 window of VALID pixels (in mm
% coords centred on the pixel) and returns dT/dx, dT/dy AT the centre plus the
% fit R^2. order=1 plane, order=2 quadratic (captures curvature). Only valid
% neighbours are used, so the irregular tissue boundary is handled correctly.
[nr, nc] = size(T);
Gx = nan(nr, nc);  Gy = nan(nr, nc);  Rsq = nan(nr, nc);
[oxg, oyg] = meshgrid(-win:win, -win:win);   % pixel offsets (col, row)
oxg = oxg(:);  oyg = oyg(:);
minpts = 3*(order >= 2) + 3;                 % 6 for quadratic, 3 for plane; need > params
for r = 1:nr
    for c = 1:nc
        if ~valid(r, c), continue; end
        rr = r + oyg;  cc = c + oxg;
        in = rr >= 1 & rr <= nr & cc >= 1 & cc <= nc;
        lin = sub2ind([nr nc], rr(in), cc(in));
        keep = valid(lin);
        lin = lin(keep);
        if numel(lin) < minpts, continue; end
        z  = T(lin);
        xm = (cc(in));  xm = (xm(keep) - c) * dx;   % mm, centred at the pixel
        ym = (rr(in));  ym = (ym(keep) - r) * dy;
        if order >= 2
            A = [ones(numel(z),1), xm, ym, xm.^2, ym.^2, xm.*ym];
        else
            A = [ones(numel(z),1), xm, ym];
        end
        coef = A \ z;                 % least-squares surface
        Gx(r, c) = coef(2);           % dT/dx at the centre (x=y=0)
        Gy(r, c) = coef(3);           % dT/dy at the centre
        resid  = z - A*coef;
        ss_tot = sum((z - mean(z)).^2);
        Rsq(r, c) = 1 - sum(resid.^2) / max(ss_tot, eps);
    end
end
end

function s = setdefault(s, field, val)
if ~isfield(s, field) || isempty(s.(field)), s.(field) = val; end
end
