function data_averaged = ensembleAverageFull(data, pacing)
    % data    - [rows x cols x frames] optical mapping array
    % pacing  - pacing stimulus vector; pass [] to auto-detect peaks from data

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

     % ---------- outlier rejection ----------
     meanAP = mean(data_matrix, 3, 'omitnan');
     rmsDev = sqrt(mean((data_matrix - meanAP).^2, 4, 'omitnan'));
     mu = mean(rmsDev, 3, 'omitnan');
     sig = std(rmsDev, 0, 3, 'omitnan');
     bad    = rmsDev > mu + 3*sig;
     bad_signals = repmat(bad,[1 1 1 size(data_matrix, 4)]);
     data_matrix(bad_signals) = NaN;  % Remove outliers by setting them to NaN
     

      % ---------- ensemble average ----------
    data_averaged = normalize_data(squeeze(mean(data_matrix, 3, 'omitnan')));

end


