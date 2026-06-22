function reentry = check_for_wf_reentry(wf, wavefronts, df, frame_rate)
% check_for_wf_reentry  Detect whether a wavefront constitutes a re-entrant circuit.
%
% Returns reentry = 1 if the wavefront shows evidence of re-entry, 0 otherwise.
%
%   wf         : tracked Wavefront object (must have lifespan and path set)
%   wavefronts : full wavefront cell array (1×num_frames)
%   df         : dominant frequency in Hz (from wf_dynamics.df)
%   frame_rate : acquisition frame rate in frames/sec
%
% Re-entry criterion: the wavefront must persist for at least 0.8 of one
% rotation cycle. cycle_length is in FRAMES = round(frame_rate / df).
% Do NOT use round(1000/df) — that gives milliseconds, not frames, and is
% only correct at exactly 1000 fps.

reentry = 0;

if nargin < 4 || isempty(frame_rate) || frame_rate <= 0
    % frame_rate required for correct frame-domain threshold — skip check
    return;
end

cycle_length = round(frame_rate / df);  % frames per rotation cycle

if wf.lifespan(1,2) < 0.8 * cycle_length
    return;
end

% --- Re-entry detection (to be implemented) ---
% A wavefront is re-entrant if a representative point (e.g. midpoint) recurs
% at approximately the same location after one cycle period.
%
% Suggested approach:
%   mid_idx   = round(wf.length(1) / 2);
%   mid_point = round(wf.location(mid_idx, :));   % [x, y]
%   num_frames  = size(wf.path, 1);
%   start_frame = round(0.8 * cycle_length);
%   for i = start_frame : num_frames
%       loc = round(wavefronts{1, wf.path(i,1)}{wf.path(i,2), 1}.location);
%       if any(loc(:,1) == mid_point(1) & loc(:,2) == mid_point(2))
%           reentry = 1;
%           return;
%       end
%   end

end
