function new_data = remove_Drift(data)
%% remove_Drift removes temporal baseline drift from optical mapping data.
%
% INPUTS
%   data     - [rows x cols x frames] optical mapping array
%
% OUTPUT
%   new_data - drift-corrected array, same size as input
%
% METHOD
%   Fits and subtracts a degree-3 polynomial from each pixel time series
%   using MATLAB's detrend().  All pixels are processed in a single
%   vectorised call — no per-pixel loop.

%% Code
% Reshape to [T x pixels] so detrend (which operates on columns) handles
% all pixel time series in one call and works on all MATLAB versions.
% permute([3,1,2]) puts time first so the subsequent reshape is a
% zero-copy view; permute([2,3,1]) reverses it at the end.
[R, C, T] = size(data);
tmp      = reshape(permute(data, [3,1,2]), T, R*C);
tmp      = detrend(tmp, 3);
new_data = permute(reshape(tmp, T, R, C), [2,3,1]);
