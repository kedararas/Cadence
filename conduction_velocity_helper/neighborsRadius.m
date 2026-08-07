function [nbr_rc, nbr_lin] = neighborsRadius(pos, sz, R)
% neighborsRadius  Flexible 2D neighborhood for any radius R.
%
%   [nbr_rc, nbr_lin] = neighborsRadius(pos, sz, R, ...)
%
% INPUTS
%   pos : [r c] or scalar linear index into array of size sz
%   sz  : [nRows nCols]
%   R   : neighborhood radius (positive integer)
%
% NAME-VALUE OPTIONS (all optional)
%   'Metric'        : 'chebyshev' (default), 'cityblock', or 'euclidean'
%   'Wrap'          : true/false (default false) — toroidal wrapping
%   'IncludeCenter' : true/false (default false)
%
% OUTPUTS
%   nbr_rc  : Kx2 [row col] list of neighbors (K depends on metric/edges)
%   nbr_lin : Kx1 linear indices into an array of size sz
%
% NOTES
%   - With Metric='chebyshev' and IncludeCenter=false, returns (2R+1)^2 - 1
%     neighbors (subject to clipping unless 'Wrap' is true).

    % arguments
    %     pos {mustBeNumeric, mustBeVector}
    %     sz  (1,2) {mustBeInteger, mustBePositive}
    %     R   (1,1) {mustBeInteger, mustBeNonnegative}
    %     varargin.Metric        (1,1) string {mustBeMember(varargin.Metric, ["chebyshev","cityblock","euclidean"])} = "chebyshev"
    %     varargin.Wrap          (1,1) logical = false
    %     varargin.IncludeCenter (1,1) logical = false
    % end

    nRows = sz(1); nCols = sz(2);

    % Normalize input to (r,c)
    if isscalar(pos)
        [r, c] = ind2sub([nRows nCols], pos);
    elseif numel(pos) == 2
        r = pos(1); c = pos(2);
    else
        error('pos must be [r c] or a scalar linear index.');
    end

    % Build offset grid
    rngv = -R:R;
    [dc, dr] = meshgrid(rngv, rngv);  % dc = column offsets, dr = row offsets
    dr = dr(:); dc = dc(:);

    % Center mask
    isCenter = (dr == 0) & (dc == 0);

    % % Metric mask
    % switch varargin.Metric
    %     case "chebyshev"
    %         inMetric = max(abs(dr), abs(dc)) <= R;
    %     case "cityblock"
    %         inMetric = (abs(dr) + abs(dc)) <= R;
    %     case "euclidean"
    %         inMetric = (dr.^2 + dc.^2) <= R^2;
    % end

    inMetric = max(abs(dr), abs(dc)) <= R;

    % Include/exclude center
    % if varargin.IncludeCenter
    %     keep = inMetric;
    % else
        keep = inMetric & ~isCenter;
    %end

    dr = dr(keep); dc = dc(keep);

    % Candidate neighbor coordinates
    cand = [r + dr, c + dc];

    % if varargin.Wrap
    %     % Toroidal wrap to 1..N
    %     cand(:,1) = mod(cand(:,1)-1, nRows) + 1;
    %     cand(:,2) = mod(cand(:,2)-1, nCols) + 1;
    % else
        % Clip to bounds
        in = cand(:,1) >= 1 & cand(:,1) <= nRows & ...
             cand(:,2) >= 1 & cand(:,2) <= nCols;
        cand = cand(in,:);
    %end

    nbr_rc  = cand;
    if isempty(cand)
        nbr_lin = zeros(0,1);
    else
        nbr_lin = sub2ind([nRows nCols], cand(:,1), cand(:,2));
    end
end
