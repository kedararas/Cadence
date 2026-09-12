function lat = compute_lat_50(stack)
            % Sub-frame activation time = frame where the OAP upstroke crosses 50% of its
            % (baseline->peak) amplitude, by linear interpolation around max dV/dt.
            % The search is CONFINED to the depolarization phase: the averaged window is
            % ~1.1x the beat cycle (ensembleAverageFull), so its second half is diastolic
            % baseline. Searching the whole trace lets max(dV/dt) latch a noise spike in
            % that tail for low-SNR pixels -> spurious late "next-beat" LATs. Pixels with
            % no real upstroke in the window are returned NaN, not a forced wrong frame.
            %
            % Guards on that confinement:
            %  * Per pixel, max(dV/dt) is taken only BEFORE the pixel's own peak inside the
            %    search window, because the upstroke rises to the peak. A one-frame step
            %    after the AP (an ensemble-average seam) or the next beat's rise can then
            %    never outrank the real upstroke.
            %  * If the population beat never relaxes below 10% inside the window, the
            %    search ends before the rise out of the population's diastolic trough.
            %  * When the max-dV/dt frame is above the 50% level, LAT is interpolated
            %    between the samples that bracket the level. When it is below (a shoulder
            %    on a low-SNR upstroke), the steepest segment is extrapolated to the level,
            %    as before; that is more stable than the literal crossing, which a small
            %    shoulder near 50% can push many frames late. The extrapolation is capped
            %    at the pixel's peak. It cannot pass it anyway while max dV/dt is taken
            %    before the peak, so the cap is only a guard.
            try
                [nr, nc, nt] = size(stack);
                M    = reshape(stack, nr*nc, nt);          % pixels x time
                dvdt = diff(M, 1, 2);                       % pixels x (nt-1); dvdt(:,j) = M(:,j+1) - M(:,j)
                base = min(M, [], 2);  peak = max(M, [], 2);
                amp  = peak - base;    thr = base + 0.5*amp;

                % --- activation window from the population beat ---
                active = amp > 0.10*max(amp);              % pixels with a clear AP
                if any(active), pop = mean(M(active,:), 1, 'omitnan');
                else,           pop = mean(M, 1, 'omitnan'); end
                pop = (pop - min(pop)) / (max(pop) - min(pop) + eps);
                [~, pk]    = max(pop);                      % AP peak
                tail       = find(pop(pk:end) < 0.10, 1);   % first return toward baseline
                if isempty(tail)
                    [~, trough] = min(pop(pk:end));         % trough frame = pk + trough - 1
                    search_end  = pk + trough - 2;          % last dV/dt sample before the rise out of it
                else
                    search_end  = pk + tail;
                end
                search_end = min(max(search_end, 1), nt-1);

                % --- max dV/dt within [1, search_end], and before each pixel's own peak ---
                [~, ppk] = max(M(:, 1:search_end+1), [], 2);
                dvs = dvdt;  dvs(:, search_end+1:end) = -Inf;
                dvs((1:nt-1) >= ppk) = -Inf;                % rise INTO the peak is dvdt(:, ppk-1)
                [slope, fmax] = max(dvs, [], 2);
                slope(isinf(slope)) = NaN;                  % peak at frame 1: no upstroke in the window

                % --- validity: real upstroke needs amplitude AND slope above noise ---
                slope_ref = median(slope(active), 'omitnan');
                good = active & (slope > 0.20*slope_ref);

                lat = nan(nr*nc, 1);                        % invalid -> NaN
                gi  = find(good);
                for ii = 1:numel(gi)
                    p = gi(ii);  v = M(p,:);  f = fmax(p);  t = thr(p);
                    k = f;
                    while k > 1 && v(k) > t,  k = k - 1;  end
                    if v(k) <= t && k < nt && v(k+1) > v(k)    % bracketed above f, or the steepest segment below 50%
                        lat(p) = min(k + (t - v(k)) / (v(k+1) - v(k)), ppk(p));
                    else
                        lat(p) = f;                            % upstroke starts at/before window edge -> max-dV/dt frame
                    end
                end
                lat = reshape(lat, nr, nc);
            catch ME
                dd = diff(stack,1,3);  [~, lat] = max(dd, [], 3);  lat = double(lat);
            end
        end
