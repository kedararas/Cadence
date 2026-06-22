function phase_singularity_data = extract_phase_singularity(phase_data)
%EXTRACT_PHASE_SINGULARITY  Locate phase singularities frame-by-frame.
%
%   phase_singularity_data = extract_phase_singularity(phase_data)
%
%   Top-level phase-singularity detector (called from the Feature Extraction
%   module). For each frame, computes the spatial phase gradient and the line
%   integral of phase around every pixel's 3x3 neighbourhood; a non-zero
%   topological charge (the phase winding +/-2*pi) marks a phase singularity --
%   the pivot of a reentrant rotor. All frame-invariant index/weight grids are
%   hoisted out of the loop so only the per-frame phase values vary.
%
%   Inputs
%     phase_data  rows x cols x num_frames real array of wrapped phase
%                 (radians), e.g. from a Hilbert-transform phase map.
%
%   Output
%     phase_singularity_data  rows x cols x num_frames array; non-zero
%                             entries flag phase-singularity pixels (sign
%                             gives chirality), zero elsewhere.

frame_len = size(phase_data, 3);
rows      = size(phase_data, 1);
cols      = size(phase_data, 2);

% Pre-allocate output and padded gradient arrays — done once, not per frame
phase_singularity_data = zeros(rows, cols, frame_len);
kbx = zeros(rows+2, cols+2);
kby = zeros(rows+2, cols+2);

%% Hoist all constant index/weight arrays outside the loop.
% These depend only on frame dimensions, not frame content — identical every frame.

% Linear index grids for the padded arrays
indX        = reshape(1:(rows+2)*(cols+2), [rows+2, cols+2]);
indInsideX  = repmat(indX(2:end-1, 2:end-1), [1, 1, 9]);

indY        = reshape(1:(rows+2)*(cols+2), [rows+2, cols+2]);
indInsideY  = repmat(indY(2:end-1, 2:end-1), [1, 1, 9]);

% 3×3 neighbor offsets for x and y gradient arrays
pb = rows + 2; % padded number of rows (stride for column-major indexing)

kLocX_kernel = [-(pb+1), -1, (pb-1); -pb, 0, pb; -(pb-1), 1, (pb+1)];
kLocX_kernel = reshape(kLocX_kernel, [1, 1, 9]);
kLocX        = repmat(kLocX_kernel, [rows, cols, 1]) + indInsideX;

kLocY_kernel = [-(pb+1), -1, (pb-1); -pb, 0, pb; -(pb-1), 1, (pb+1)];
kLocY_kernel = reshape(kLocY_kernel, [1, 1, 9]);
kLocY        = repmat(kLocY_kernel, [rows, cols, 1]) + indInsideY;

% Sobel-like quadrature weights for topological charge computation
COx = reshape([-0.5, 0, 0.5; -1, 0, 1; -0.5, 0, 0.5], [1, 1, 9]);
COx = repmat(COx, [rows, cols, 1]);

COy = reshape([0.5, 1, 0.5; 0, 0, 0; -0.5, -1, -0.5], [1, 1, 9]);
COy = repmat(COy, [rows, cols, 1]);

%% Main loop — only frame-varying work remains inside
for i = 1:frame_len
    phase = phase_data(:,:,i);

    % x-axis partial derivative with wrap-around boundary
    kx = [phase(:,2:end) - phase(:,1:end-1), phase(:,end)];
    kbx(2:end-1, 2:end-1) = kx;
    kbx = mod(kbx + pi, 2*pi) - pi;   % wrap to [-pi, pi] in one pass

    % y-axis partial derivative with wrap-around boundary
    ky = [phase(2:end,:) - phase(1:end-1,:); phase(end,:)];
    kby(2:end-1, 2:end-1) = ky;
    kby = mod(kby + pi, 2*pi) - pi;

    % Topological charge via discrete line integral
    nt = sum(kby(kLocY) .* COx, 3) + sum(kbx(kLocX) .* COy, 3);

    % Normalize by the dominant extremum
    nt_max = max(nt(:));
    nt_min = min(nt(:));
    if abs(nt_max) > abs(nt_min)
        nt = nt / nt_max;
    else
        nt = nt / nt_min;
    end

    phase_singularity_data(:,:,i) = nt;
end

end

