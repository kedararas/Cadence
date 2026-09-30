function data_averaged = ensembleAverageFull(data, pacing)
    % data    - [rows x cols x frames] optical mapping array
    % pacing  - pacing stimulus vector; pass [] to auto-detect peaks from data
    %
    % Per pixel, beats are baseline-corrected, outlier beats are dropped by a
    % robust median/MAD rule on their RMS deviation from the mean beat (>= 4
    % beats; see below), and the rest are averaged and normalized to [0, 1].

    [X, Y, T] = size(data);

    if isempty(pacing)
        pacing = auto_detect_peaks(data);
    end

    % Support both binary onset vectors (from auto_detect_peaks, which places
    % a 1 at each beat onset) and smooth analog pacing traces (findpeaks).
    % For a binary vector, find() is exact and avoids findpeaks edge cases
    % (e.g. a value of 1 adjacent to another 1 is not a strict local max).
    if all(pacing == 0 | pacing == 1)
        locs = find(pacing);
    else
        [~, locs] = findpeaks(pacing);
    end
    % Use median IBI for num_frames so recordings where pacing does not
    % span the full duration (e.g. only 2 beats, starting late) get the
    % correct beat length rather than T/n_beats.
    if length(locs) >= 2
        num_frames = round(median(diff(locs)));
    else
        num_frames = round(T / length(locs));
    end
    pre_window   = round(num_frames * 0.1);
    total_frames = pre_window + num_frames;

    n_beats     = length(locs);
    data_matrix = NaN(X, Y, n_beats, total_frames);

    % ---------- extract windows ----------
    % Beats at the start/end of the recording may not have a full pre-window
    % or a full post-onset window.  Rather than discarding them entirely,
    % copy whatever frames ARE available into the correct position inside the
    % NaN-initialised window.  The omitnan average later handles missing edges.
    valid_beats = false(1, n_beats);
    for i = 1:n_beats
        i_start = locs(i) - pre_window;    % may be < 1
        i_end   = locs(i) + num_frames - 1;% may be > T

        % Clamp to recording bounds
        src_start = max(i_start, 1);
        src_end   = min(i_end,   T);

        if src_end < src_start
            continue;                       % beat onset is outside recording
        end

        % Destination positions inside the total_frames window
        dst_start = src_start - i_start + 1;   % 1-based offset
        dst_end   = dst_start + (src_end - src_start);

        data_matrix(:, :, i, dst_start:dst_end) = data(:, :, src_start:src_end);
        valid_beats(i) = true;
    end
    data_matrix = data_matrix(:, :, valid_beats, :);

    if size(data_matrix, 3) == 0
        warning('ensembleAverageFull:noValidBeats', ...
            'No beat windows fell within the recording bounds. Returning zeros.');
        data_averaged = zeros(X, Y, total_frames);
        return;
    end

    % ---------- per-beat baseline removal ----------
    baseline = mean(data_matrix(:,:,:,1:pre_window), 4, 'omitnan');  % [X Y beats]
    data_matrix = data_matrix - reshape(baseline, [X Y size(data_matrix,3) 1]);

    % ---------- outlier rejection (per pixel, robust) ----------
    % Each beat's RMS deviation from the mean beat is judged against the
    % MEDIAN of those deviations across beats, with a robust spread
    %     rsd = max(1.4826 * MAD, 0.10 * median)
    % and the beat is dropped at that pixel when rmsDev > median + 3 * rsd.
    %
    % Why not mean + 3*SD (the previous rule): the outlier is part of the mean
    % and SD it is judged against. With the N-1 SD the largest attainable z is
    % (n-1)/sqrt(n) -- 2.85 at n = 10 -- so the rule could never reject a beat
    % in a recording of 10 or fewer beats, however bad the beat was.
    %
    % The 10%-of-median floor on rsd stops near-identical beats (MAD -> 0)
    % from being rejected over trivial differences. In synthetic tests (250-
    % frame beats, noise 0.05-0.3 of amplitude) it gave 0% rejection of clean
    % beats at every n, caught a half-beat 0.5-amplitude artifact in >=91% of
    % pixels for n >= 5, and did not reject 15%-amplitude alternans beats for
    % n >= 5. A 0.25 floor protected alternans further but missed that
    % artifact entirely at noise 0.3.
    %
    % Needs >= 4 beats: with 3 the median and MAD are set by single beats, and
    % the minority beat of an alternating pair was rejected in ~30% of pixels
    % at high SNR. At most floor(n/2) beats can exceed the median, so a pixel
    % never loses all its beats.
    MIN_BEATS_REJECT = 4;
    if size(data_matrix, 3) >= MIN_BEATS_REJECT
        meanAP = mean(data_matrix, 3, 'omitnan');
        rmsDev = sqrt(mean((data_matrix - meanAP).^2, 4, 'omitnan'));   % [X Y beats]
        med    = median(rmsDev, 3, 'omitnan');
        rsd    = max(1.4826 * median(abs(rmsDev - med), 3, 'omitnan'), 0.10 * med);
        bad    = rmsDev > med + 3*rsd;                  % NaN pixels compare false
        data_matrix(repmat(bad, [1 1 1 size(data_matrix, 4)])) = NaN;
    end


      % ---------- ensemble average ----------
    data_averaged = normalize_data(squeeze(mean(data_matrix, 3, 'omitnan')));

end


