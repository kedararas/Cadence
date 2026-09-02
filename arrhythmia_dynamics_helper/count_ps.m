function [ps, ps_count] = count_ps(ps_data, wavefronts, pixel_size_mm, ...
                                   ps_charge_threshold, min_opposite_sep_px)
% count_ps  Detect, deduplicate, and annotate phase singularities per frame.
%
% Wavefront proximity is treated as METADATA, not a hard filter.
% ps{1,i} contains all cleaned PS with wavefront association columns appended.
% ps{2,i} contains all cleaned PS (no wavefront info).
% ps{3,i} contains all raw PS before deduplication.
%
% Downstream code can filter by col 5 > 0 to isolate WF-associated PS,
% or use all rows and treat col 5 as descriptive metadata.
%
% Optional argument pixel_size_mm (mm/pixel) enables physically-grounded
% thresholds. If omitted, thresholds are derived from frame geometry.
%
% Optional argument ps_charge_threshold (default 3.5) is the absolute
% topological-charge cut used to accept a pixel as a PS (see below).
%
% Optional argument min_opposite_sep_px (default 3) is the closest two
% OPPOSITE-charge PS may sit and still both be kept. It is deliberately in
% PIXELS and does NOT scale with pixel_size_mm: it suppresses a discretisation
% artefact of the charge operator, whose footprint is set by the 3x3 stencil,
% not by tissue geometry (see the deduplication block below).

num_frames = size(ps_data, 3);
num_rows   = size(ps_data, 1);
num_cols   = size(ps_data, 2);
frame_diag = sqrt(num_rows^2 + num_cols^2);

% Absolute topological-charge threshold for accepting a pixel as a PS.
% ps_data now carries the RAW winding charge (see extract_phase_singularity),
% where a full +/-2*pi singularity reads +/-6.283 -- verified by running the
% operator on an analytic spiral, atan2(y-y0, x-x0), which returns exactly
% +2*pi at the core -- and background noise sits below ~2. (An earlier note
% here claimed ~8.9; that was measured before the per-frame normalisation was
% removed and does not describe the raw charge.) A single core spreads charge
% over its neighbours, so cluster maxima can exceed one winding. This is a
% fixed-reference test (unlike the old
% per-frame-normalized round()==+/-1), so a given rotor core is detected
% consistently instead of blinking in/out as neighbouring PS wax and wane.
% Default 3.5 catches the falloff skirt of real cores — improving detection
% continuity — while staying above the noise floor: measured on rat VF, dropping
% 4.4 -> 3.5 added only ~0.09 genuinely-new PS/frame and left PS density flat.
% Tuning: raise toward 4.4-5 for fewer/cleaner PS; do NOT go below ~3.2, where
% the low-charge noise mass climbs steeply. VALIDATE against visual scoring
% before trusting a new value across recordings/species.
if nargin < 4 || isempty(ps_charge_threshold)
    ps_charge_threshold = 3.5;
end

% Closest two OPPOSITE-charge PS may sit and both be kept. See the
% deduplication block for why this is a pixel constant rather than a distance
% in tissue.
if nargin < 5 || isempty(min_opposite_sep_px)
    min_opposite_sep_px = 3;
end

if nargin >= 3 && ~isempty(pixel_size_mm) && pixel_size_mm > 0
    % Physical thresholds — grounded in tissue scale
    % PS detection clusters span ~1-2 mm; tip association within ~0.5 mm
    dup_threshold      = round(2.0 / pixel_size_mm);  % merge clusters within 2 mm
    wf_assoc_threshold = round(0.5 / pixel_size_mm);  % PS within 0.5 mm of wavefront tip
    fprintf(['count_ps: pixel_size=%.3f mm -> dup=%d px, wf_assoc=%d px, ' ...
             'charge_thr=%.2f, opp_sep=%d px\n'], ...
        pixel_size_mm, dup_threshold, wf_assoc_threshold, ps_charge_threshold, ...
        min_opposite_sep_px);
else
    % Geometry-derived fallbacks when pixel size is unknown
    % Deduplication: ~3% of frame diagonal — merges local detection clusters
    % WF association: ~1% of frame diagonal — PS must be close to wavefront tip
    dup_threshold      = max(5,  round(0.03 * frame_diag));
    wf_assoc_threshold = max(3,  round(0.01 * frame_diag));
    fprintf(['count_ps: no pixel_size supplied -> dup=%d px, wf_assoc=%d px, ' ...
             'charge_thr=%.2f, opp_sep=%d px\n'], ...
        dup_threshold, wf_assoc_threshold, ps_charge_threshold, min_opposite_sep_px);
end

% Pre-allocate outputs — 3 rows each
ps       = cell(3, num_frames);
ps_count = zeros(3, num_frames);

for i = 1:num_frames
    nt_frame = squeeze(ps_data(:,:,i));

    % Detect positive and negative topological charges by absolute winding:
    % a genuine PS clears +/-ps_charge_threshold regardless of the other PS in
    % this frame (no per-frame normalization -> no blinking).
    [p_rows, p_cols] = find(nt_frame >=  ps_charge_threshold);
    [n_rows, n_cols] = find(nt_frame <= -ps_charge_threshold);
    ps_points = [p_cols, p_rows,  ones(numel(p_rows),1); ...
                 n_cols, n_rows, -ones(numel(n_rows),1)];

    if isempty(ps_points)
        ps{1,i} = [];
        ps{2,i} = [];
        ps{3,i} = [];
        continue;
    end

    ps{3,i} = ps_points;

    %% Deduplicate, strongest first, with two different radii by charge sign.
    %
    % SAME sign, within dup_threshold -> duplicate. Deduplication exists to
    % merge the several adjacent detections thrown off by ONE core, which
    % necessarily share that core's sign; it is not a noise filter (the charge
    % threshold above is). A sign-blind radius instead deleted real objects:
    % because ps_points lists every positive before every negative, it always
    % resolved mixed pairs in the positives' favour -- measured on rat VF,
    % positives survived 100% and negatives 12.7%, so whether a rotor could be
    % tracked depended on which way it span. choose_ps_cand enforces charge
    % conservation when linking frames, so those censored negatives broke the
    % tracks of real negative-chirality rotors (core held in 38-56% of frames,
    % vs 59-77% once sign-aware).
    %
    % OPPOSITE sign, within min_opposite_sep_px -> unresolvable, keep the
    % stronger. Making the radius sign-aware exposed an artefact the old
    % sign-blind rule had been masking: the discrete winding operator responds
    % to the phase-WRAP LINE as well as to its endpoints, emitting a dipole --
    % adjacent lobes of opposite sign, both well over threshold (measured
    % +5.5 next to -5.2 two pixels apart). Those are one feature, not two
    % cores; the phase cannot wind twice inside two pixels. On rat VF they were
    % 45% of all detections, they are NOT removed by the downstream lifespan
    % cut (they join >=60-frame tracks at the same 10% rate as everything
    % else, because the wrap line itself persists), and the measured
    % nearest-opposite-neighbour distribution separates cleanly: a spike below
    % 2.5 px on the pixel lattice, then flat from 3.5 px out.
    %
    % That separation is why min_opposite_sep_px is a PIXEL constant and is
    % NOT derived from pixel_size_mm. The lobe spacing is set by the operator's
    % 3x3 stencil, so it is fixed in pixels whatever a pixel measures; a
    % physical radius would drift with the FOV for no reason connected to the
    % artefact. On this rat prep (~0.067 mm/px) a 1 mm rule would have cut at
    % 15 px and taken out most of the genuine distribution.
    %
    % Both tests need a defensible winner, so candidates are visited in
    % descending |charge| and the surviving detection in each cluster is its
    % strongest member -- the best available estimate of the core. (The old
    % loop kept whichever came first in raster order.)
    num_ps    = size(ps_points, 1);
    distances = pdist2(ps_points(:,1:2), ps_points(:,1:2));
    same_sign = ps_points(:,3) == ps_points(:,3)';   % [num_ps x num_ps] logical
    ps_mag    = abs(nt_frame(sub2ind([num_rows num_cols], ...
                                     ps_points(:,2), ps_points(:,1))));
    [~, order] = sort(ps_mag, 'descend');

    keep     = false(num_ps, 1);
    kept_idx = zeros(num_ps, 1);
    n_kept   = 0;
    for ii = 1:num_ps
        x = order(ii);
        if n_kept > 0
            k = kept_idx(1:n_kept);
            too_close = ( same_sign(x,k)' & (distances(x,k)' < dup_threshold)) | ...
                        (~same_sign(x,k)' & (distances(x,k)' < min_opposite_sep_px));
            if any(too_close), continue; end
        end
        n_kept           = n_kept + 1;
        kept_idx(n_kept) = x;
        keep(x)          = true;
    end
    cleaned_ps_points = ps_points(keep, :);

    %% Tag each PS to its nearest wavefront (metadata only — not a filter)
    wf_tagged = tag_ps_to_wf(cleaned_ps_points, wavefronts{2,i}, wf_assoc_threshold);

    ps{1,i} = wf_tagged;
    ps{2,i} = cleaned_ps_points;

    ps_count(1,i) = sum(wf_tagged(:,5) > 0);       % PS with a wavefront association
    ps_count(2,i) = size(cleaned_ps_points, 1);     % all cleaned PS
    ps_count(3,i) = size(ps_points, 1);             % raw PS before dedup
end
