function [B lims] = rescale_sat(A, lims, out_lims)
%RESCALE_SAT  Linearly rescale values in an array, saturating out-of-range values.
%
% Examples:
%   B = rescale_sat(A)
%   B = rescale_sat(A, lims)
%   B = rescale_sat(A, lims, out_lims)
%   [B lims] = rescale_sat(A)
%
% Linearly rescales values in an array, saturating values outside limits.
%
% NOTE: Renamed from 'rescale' to avoid shadowing MATLAB's built-in rescale,
% which has a DIFFERENT signature (rescale(A, lo, hi) sets output bounds only)
% and version-dependent handling of constant/empty frames. This custom version
% takes an explicit input saturation range and guards the divide-by-zero case
% by returning zeros, which several display paths (real2rgb) rely on.
%
% IN:
%   A - Input array of any size and class.
%   lims - 1x2 array of saturation limits to be used on A. Default:
%          [min(A(:)) max(A(:))].
%   out_lims - 1x2 array of output limits the values in lims are to be
%              rescaled to. Default: [0 1].
%
% OUT:
%   B - size(A) double array.
%   lims - 1x2 array of saturation limits used on A. Equal to the input
%          lims, if given.

% Copyright: Oliver Woodford, 2009 - 2011
% Source: MATLAB Central File Exchange (Oliver Woodford). Renamed from 'rescale'
% and lightly modified for CADENCE. License: BSD-2-Clause (see THIRD_PARTY_LICENSES.md).

if nargin < 3
    out_lims = [0 1];
end
if nargin < 2 || isempty(lims)
    M = isfinite(A);
    if ~any(reshape(M, numel(M), 1))
        % All NaNs, Infs or -Infs
        B = double(A > 0);
        lims = [0 1];
    else
        lims = [min(A(M)) max(A(M))];
        B = normalize(A, lims, out_lims);
        B = min(max(B, out_lims(1)), out_lims(2));
    end
    clear M
else
    B = normalize(A, lims, out_lims);
    B = min(max(B, out_lims(1)), out_lims(2));
end
return

function B = normalize(A, lims, out_lims)
if lims(2) == lims(1) || out_lims(1) == out_lims(2)
    B = zeros(size(A));
else
    B = double(A);
    if lims(1)
        B = B - lims(1);
    end
    v = (out_lims(2) - out_lims(1)) / (lims(2) - lims(1));
    if v ~= 1
        B = B * v;
    end
end
if out_lims(1)
    B = B + out_lims(1);
end
return
