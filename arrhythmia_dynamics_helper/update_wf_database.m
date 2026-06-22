function wf_dynamics = update_wf_database(wavefront_dynamics, wf, wavefronts, threshold)
% WavefrontDynamics is a handle class — this is a reference, not a copy.
wf_dynamics = wavefront_dynamics;

% Cache repeated expressions computed from wf fields
wf_bd1          = wf.birthday(1,1);
wf_bd2          = wf.birthday(1,2);
wf_end_length   = wf.length(end,1);
wf_lifespan     = wf.lifespan(1,2);
wf_median_len   = round(median(wf.length));
thr7            = threshold(1,7);
thr8            = threshold(1,8);

% Add entry for wf size & duration
% end+1 on an empty matrix correctly appends as row 1 — no isempty guard needed
wf_dynamics.wf_size_duration(end+1,:) = [wf_bd1, wf_bd2, wf_median_len, wf_end_length, wf_lifespan];

% Cache current size_duration row count used as record ID throughout
wf_sd_index = size(wf_dynamics.wf_size_duration, 1);
base_row    = [wf_sd_index, wf_bd1, wf_bd2, wf_end_length, wf_lifespan];

% Add entry for wf path (long-lived, large wavefronts only)
if wf_lifespan >= thr8 && wf_median_len >= 6*thr7
    n = size(wf_dynamics.wf_path, 1);
    wf_dynamics.wf_path{n+1, 1} = [wf_sd_index, wf_bd1, wf_bd2, wf_median_len, wf_lifespan];
    wf_dynamics.wf_path{n+1, 2} = wf.path;
end

% Add cause-of-death entry using switch to avoid repeated strcmp calls
switch wf.cause_of_death
    case 'FRAGMENTED'
        wf_dynamics.wf_fractionations(end+1,:) = base_row;
    case 'MERGED'
        wf_dynamics.wf_collisions(end+1,:) = base_row;
    case 'EXPIRED'
        wf_dynamics.wf_blocks(end+1,:) = base_row;
end

% Add breakthrough entry: wf has no parent and all locations are away from
% field boundaries (defined as within thr7 of any edge)
if isempty(wf.parent) && ...
        ~any(wf.location(:,1) > 9*thr7) && ~any(wf.location(:,2) > 9*thr7) && ...
        ~any(wf.location(:,1) < thr7)   && ~any(wf.location(:,2) < thr7)
    wf_dynamics.wf_breakthroughs(end+1,:) = base_row;
end

% Check repeatability and uniqueness
[repeat_index, f_dist, e_dist] = check_for_wf_repeatability(wf_dynamics, wf, wavefronts, threshold);

if repeat_index > 0
    % Repeated wavefront pattern
    wf_dynamics.wf_repeatability(end+1,:) = [wf_sd_index, repeat_index, wf_bd1, wf_bd2, wf_median_len, wf_lifespan, f_dist, e_dist];
elseif wf_lifespan >= thr8
    % Unique long-lived wavefront
    wf_dynamics.wf_multiplicity(end+1,:) = [wf_sd_index, wf_bd1, wf_bd2, wf_median_len, wf_lifespan];
end

% Check for reentry
reentry_index = check_for_wf_reentry(wf, wavefronts, wf_dynamics.df, wf_dynamics.frame_rate);
if reentry_index == 1
    wf_dynamics.wf_reentry(end+1,:) = [wf_sd_index, wf_bd1, wf_bd2, wf_median_len, wf_lifespan];
end
