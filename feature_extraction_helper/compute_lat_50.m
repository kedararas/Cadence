function lat = compute_lat_50(stack)
            % Sub-frame activation time = frame where the OAP upstroke crosses 50% of its
            % (baseline->peak) amplitude, by linear interpolation around max dV/dt.
            % The search is CONFINED to the depolarization phase: the averaged window is
            % ~1.1x the beat cycle (ensembleAverageFull), so its second half is diastolic
            % baseline. Searching the whole trace lets max(dV/dt) latch a noise spike in
            % that tail for low-SNR pixels -> spurious late "next-beat" LATs. Pixels with
            % no real upstroke in the window are returned NaN, not a forced wrong frame.
            try
                [nr, nc, nt] = size(stack);
                M    = reshape(stack, nr*nc, nt);          % pixels x time
                dvdt = diff(M, 1, 2);                       % pixels x (nt-1)
                base = min(M, [], 2);  peak = max(M, [], 2);
                amp  = peak - base;    thr = base + 0.5*amp;

                % --- activation window from the population beat ---
                active = amp > 0.10*max(amp);              % pixels with a clear AP
                if any(active), pop = mean(M(active,:), 1, 'omitnan');
                else,           pop = mean(M, 1, 'omitnan'); end
                pop = (pop - min(pop)) / (max(pop) - min(pop) + eps);
                [~, pk]    = max(pop);                      % AP peak
                tail       = find(pop(pk:end) < 0.10, 1);   % first return toward baseline
                if isempty(tail), search_end = nt - 1;
                else,             search_end = min(nt-1, pk + tail); end

                % --- max dV/dt within [1, search_end] only ---
                dvs = dvdt;  dvs(:, search_end+1:end) = -Inf;
                [slope, fmax] = max(dvs, [], 2);

                % --- validity: real upstroke needs amplitude AND slope above noise ---
                slope_ref = median(slope(active), 'omitnan');
                good = active & (slope > 0.20*slope_ref);

                lat = nan(nr*nc, 1);                        % invalid -> NaN
                gi  = find(good);
                for ii = 1:numel(gi)
                    p = gi(ii);  v = M(p,:);  f = fmax(p);  t = thr(p);
                    k = f;
                    while k > 1 && v(k) > t,  k = k - 1;  end
                    if v(k) <= t && k < nt && v(k+1) > v(k)    % threshold bracketed: v(k) <= t <= v(k+1)
                        lat(p) = k + (t - v(k)) / (v(k+1) - v(k));   % frac in [0,1] -> lat in [k, k+1], always >= 1
                    else
                        lat(p) = f;                            % upstroke starts at/before window edge -> max-dV/dt frame
                    end
                end
                lat = reshape(lat, nr, nc);
            catch ME
                dd = diff(stack,1,3);  [~, lat] = max(dd, [], 3);  lat = double(lat);
            end
        end