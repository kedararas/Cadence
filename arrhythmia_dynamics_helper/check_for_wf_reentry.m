function reentry = check_for_wf_reentry(wf, wavefronts, df, frame_rate, pixel_size_mm)
% check_for_wf_reentry  Detect whether a tracked wavefront constitutes functional
% reentry — i.e. it circulates and returns to tissue it occupied ~one rotation
% earlier.
%
% Returns reentry = 1 if the wavefront DEPARTS from its birth region and RETURNS
% near it within the 0.8-1.5 rotation window, 0 otherwise.
%
%   wf            : tracked Wavefront object (lifespan and path must be set)
%   wavefronts    : full 1 x num_frames wavefront cell array
%   df            : rotor frequency in Hz. Pass the SAME estimate the PS pipeline
%                   uses (estimate_rotor_frequency, p75 of the DF map) so the
%                   rotation period is consistent across the two analyses — NOT
%                   the plain median, which a bimodal map drags toward slow
%                   far-field tissue, lengthening the cycle.
%   frame_rate    : acquisition frame rate (frames/sec).
%   pixel_size_mm : (optional) mm/pixel. If given, the departure/return
%                   tolerances are physical (~2 mm / ~1 mm); otherwise they scale
%                   to the wavefront's own arc length so the test still works
%                   without a known pixel size.
%
% cycle_length is in FRAMES = round(frame_rate / df). Do NOT use round(1000/df)
% — that gives milliseconds, not frames, and is only correct at exactly 1000 fps.
%
% Why departure AND return: a wavefront in reentry sweeps around a core and
% comes back, so its centroid must leave the birth region (this rejects a
% wavefront that merely lingers in place) and then return near it about one
% rotation later. This is the lightweight criterion; a tip winding-angle test
% (accumulated rotation >= 2*pi about the curve centroid) is the rigorous
% upgrade if false positives appear. VALIDATE against known reentrant episodes
% before trusting the wf_reentry counts.

reentry = 0;

if nargin < 4 || isempty(frame_rate) || frame_rate <= 0 || isempty(df) || df <= 0
    % frame_rate and df required for the frame-domain rotation period.
    return;
end
if nargin < 5
    pixel_size_mm = [];
end
if isempty(wf.path) || size(wf.path, 1) < 3
    % Need a few tracked frames to describe a trajectory.
    return;
end

cycle_length = round(frame_rate / df);   % frames per rotation cycle
if cycle_length < 1
    return;
end
if wf.lifespan(1,2) < 0.8 * cycle_length
    % Too short-lived to close a reentrant loop.
    return;
end

% --- Centroid trajectory along the wavefront's own tracked path ---
% Each path row is [frame, index]; the per-frame curve is retrieved from the
% wavefront cell array and reduced to its centroid.
num_links = size(wf.path, 1);
centroid  = nan(num_links, 2);
for p = 1:num_links
    seg = wavefronts{1, wf.path(p,1)}{wf.path(p,2), 1}.location;  % L x 2 [x, y]
    centroid(p,:) = mean(seg, 1);
end
ref  = centroid(1,:);                                       % birth-region centroid
dist = hypot(centroid(:,1) - ref(1), centroid(:,2) - ref(2)); % excursion from birth

% --- Departure / return tolerances ---
if ~isempty(pixel_size_mm) && pixel_size_mm > 0
    depart_min = 2.0 / pixel_size_mm;   % must leave by ~2 mm (a genuine excursion)
    return_tol = 1.0 / pixel_size_mm;   % must come back within ~1 mm of birth
else
    scale = median(wf.length);          % ~wavefront arc length in px (point count)
    if ~isfinite(scale) || scale <= 0
        scale = 6;                       % minimum trackable segment length
    end
    depart_min = 0.5  * scale;
    return_tol = 0.25 * scale;
end

% --- Excursion-and-return over the ~one-rotation window ---
win_hi = min(num_links, round(1.5 * cycle_length));
ret_lo = max(2, round(0.8 * cycle_length));
ret_lo = min(ret_lo, win_hi);

if max(dist(1:win_hi)) >= depart_min && min(dist(ret_lo:win_hi)) <= return_tol
    reentry = 1;
end

end
