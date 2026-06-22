function [Vmag, Vx, Vy, Gx, Gy, meta] = conduction_velocity(T, dx, dy, opts)
% CONDUCTION_VELOCITY  Compute conduction velocity from an activation-time map.
%
% Inputs
%   T    : activation time map (2D, seconds). NaNs allowed for invalid pixels
%   dx   : spatial step in x-direction (e.g., mm per pixel)
%   dy   : spatial step in y-direction (e.g., mm per pixel)
%   opts : struct with optional fields:
%          .smooth_sigma_pix  : Gaussian sigma in pixels (default 0 = none)
%          .grad_min          : minimum |grad T| to trust (ms/cm) (default auto)
%          .speed_max         : cap CV magnitude (mm/ms) (default 2.50)
%          .use_imgaussfilt   : true to use imgaussfilt if available (default true)
%
% Outputs
%   Vmag : CV magnitude (mm/ms)
%   Vx   : CV x-component (mm/ms)
%   Vy   : CV y-component (mm/ms)
%   Gx   : dT/dx (ms/mm)
%   Gy   : dT/dy (ms/mm)
%   meta : struct with masks & settings (valid mask, thresholds, etc.)
%
% Relationship used:
%   speed  = 1 / ||∇T||,
%   vector = -∇T / ||∇T||^2
%   (direction points from early->late; units consistent if T in s and x,y in mm)

if nargin < 4, opts = struct(); end
opts = setdefault(opts, 'smooth_sigma_pix', 2);
opts = setdefault(opts, 'grad_min', []);
opts = setdefault(opts, 'speed_max', 2.50);        % adjust to your prep/species
opts = setdefault(opts, 'use_imgaussfilt', true);

T = double(T);
validT = isfinite(T);

% --- Optional smoothing (on T, not on gradients) ---
if opts.smooth_sigma_pix > 0
    if opts.use_imgaussfilt && exist('imgaussfilt','file')
        Ts = imgaussfilt(T, opts.smooth_sigma_pix, 'FilterDomain', 'spatial');
    else
        % fallback: separable Gaussian via conv2
        sig = opts.smooth_sigma_pix;
        ksz = max(3, 2*ceil(3*sig)+1);
        g = exp(-((-(ksz-1)/2:(ksz-1)/2).^2)/(2*sig^2));
        g = g / sum(g);
        Ts = conv2(conv2(T, g, 'same'), g', 'same');
    end
else
    Ts = T;
end

% Preserve NaN regions
Ts(~validT) = NaN;

% --- Gradients (MATLAB: [dT/dy, dT/dx] = gradient(T, dy, dx)) ---
%[Gy, Gx] = gradient(Ts, dy, dx);  % Gy=dT/dy (s/cm), Gx=dT/dx (s/cm)
[Gx, Gy] = gradient(Ts, dy, dx);  % Gy=dT/dy (s/cm), Gx=dT/dx (s/cm)

% --- Gradient magnitude and validity ---
gradMag = hypot(Gx, Gy);          % ||∇T|| (s/cm)
% Auto threshold if not provided: small fraction of robust central tendency
if isempty(opts.grad_min)
    gm = gradMag(validT);
    gm = gm(isfinite(gm));
    if isempty(gm)
        opts.grad_min = 0;
    else
        opts.grad_min = max(1e-6, 0.05*median(gm));  % 5% of median
    end
end

good = validT & isfinite(gradMag) & (gradMag >= 0.5*opts.grad_min);

% --- Conduction velocity (cm/s) ---
Vmag = nan(size(T));
Vx   = nan(size(T));
Vy   = nan(size(T));

Vmag(good) = 1 ./ gradMag(good);
% Vector form: v = -∇T / ||∇T||^2
%Vx(good)   = -Gx(good) ./ (gradMag(good).^2);
Vx(good)   = Gx(good) ./ (gradMag(good).^2);
Vy(good)   = Gy(good) ./ (gradMag(good).^2);
%Vy(good)   = -Gy(good) ./ (gradMag(good).^2);

% Cap unphysiological speeds (optional)
if ~isempty(opts.speed_max) && isfinite(opts.speed_max)
    tooFast = Vmag > opts.speed_max;
    Vmag(tooFast) = opts.speed_max;
    scale = opts.speed_max ./ hypot(Vx, Vy);
    scale(~isfinite(scale)) = 1;
    Vx(tooFast) = Vx(tooFast) .* scale(tooFast);
    Vy(tooFast) = Vy(tooFast) .* scale(tooFast);
end

meta = struct();
meta.validT     = validT;
meta.good       = good;
meta.grad_min   = opts.grad_min;
meta.speed_max  = opts.speed_max;
meta.dx         = dx;
meta.dy         = dy;
meta.smooth_sigma_pix = opts.smooth_sigma_pix;

end

function s = setdefault(s, field, val)
if ~isfield(s, field) || isempty(s.(field)), s.(field) = val; end
end
