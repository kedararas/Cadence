function normData = normalize_data(data, pct)
%% Normalizes CMOS data to [0,1] using robust percentile bounds.
%
% INPUTS
% data = cmos data [rows x cols x frames]
% pct  = [lo hi] percentiles used as the normalization bounds (default [1 99]).
%        Pass [0 100] to reproduce the original min-max behaviour.
%
% OUTPUT
% normData = normalized data, same size as input, clipped to [0,1]
%
% METHOD
% For each pixel the lower/upper bounds are the lo/hi percentiles of that
% pixel's time course rather than its raw min/max.  A single outlier sample
% -- a startup spike, a cosmic ray, a filter edge transient -- can therefore
% no longer define the range and compress the real signal into a sliver of
% [0,1].  Each pixel is mapped (x - lo) / (hi - lo) and clipped to [0,1], so
% samples beyond the percentile band saturate at 0 or 1.  With pct = [0 100]
% this reduces exactly to the previous min-max normalization.

%% Code
if nargin < 2 || isempty(pct), pct = [1 99]; end

% Per-pixel percentile bounds over time (dim 3): result is [R x C x 2].
b      = prctile(data, pct, 3);
lo     = b(:,:,1);                       % [R x C] lower bound
drange = b(:,:,2) - lo;                  % [R x C] robust range
drange(drange == 0) = 1;                 % constant pixels -> numerator 0 -> stays 0

% Implicit broadcasting expands lo/drange across T without a full copy.
normData = (data - lo) ./ drange;
% Clip outliers to [0,1] WITHOUT destroying NaN.  Do not use min(max(.,0),1):
% MATLAB's max(NaN,0) returns 0, so that would convert masked-out (NaN)
% background pixels into 0.  Comparisons against NaN are false, so indexing
% leaves NaN untouched.
normData(normData < 0) = 0;
normData(normData > 1) = 1;

end
