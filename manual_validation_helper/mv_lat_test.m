function mv_lat_test()
%MV_LAT_TEST  compute_lat_50 against analytic activation times.
%
%   mv_lat_test()
%
%   Regression test, not a study instrument: synthetic stacks whose 50%
%   upstroke crossing is known exactly, one per failure mode found in the
%   2026-09-11 manual validation.
%
%     A. NEXT-BEAT LATCH.  The population beat never relaxes below 10% inside
%        the window and the next beat's (steeper) upstroke starts before the
%        window ends.  The search must stop at the diastolic trough; before
%        the fix max(dV/dt) latched the next beat and LAT landed at the end
%        of the window (rabbit 4-18: 1.4% of pixels).
%     B. MAX dV/dt BELOW THE 50% LEVEL.  A fast foot followed by a slower
%        upper upstroke.  By design (decided 2026-09-11) the steepest segment
%        is extrapolated to the 50% level rather than walking to the literal
%        crossing, which low-SNR shoulders push many frames late.  This locks
%        that behaviour; it is not a bug detector.
%     C. ORDINARY UPSTROKE.  A smooth sigmoid that relaxes fully: sub-frame
%        accuracy must hold (guards against the fixes costing precision).
%     D. ONE-FRAME SEAM AFTER THE AP.  The mechanism actually seen in rabbit
%        4-18: the population relaxes below 10% only on the last frame, so the
%        search spans the window, and a one-frame step near its end (an
%        ensemble-average seam) is steeper than the pixel's own upstroke.
%        max(dV/dt) must be taken before the pixel's own peak.
%
%   Piecewise-linear upstrokes (A, B) are sampled exactly by linear
%   interpolation, so those checks use a 1e-9 frame tolerance.

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..', 'feature_extraction_helper'));

    fprintf('mv_lat_test: compute_lat_50 against analytic activation times\n');
    nt = 165;  k = 1:nt;
    [R, C] = ndgrid(1:12, 1:12);
    ta = 30.25 + 0.37 * (C - 1) + 0.11 * (R - 1);        % sub-frame onsets, 30.25-35.5

    % ---- A. next-beat latch -------------------------------------------------
    % Upstroke 0 -> 1 over 10 frames from ta, a 2-frame plateau at exactly 1
    % (so the sampled peak is 1 and the 50% level is exactly 0.5), decay to a
    % 0.15 floor (never below 10%), then the next beat rises at 0.25/frame
    % from frame 150.
    X = zeros([size(ta) nt]);
    for i = 1:numel(ta)
        [r, c] = ind2sub(size(ta), i);
        v = min(max((k - ta(i)) / 10, 0), 1);            % upstroke, 0.1/frame
        dec = k >= ta(i) + 12;
        v(dec) = 0.15 + 0.85 * exp(-(k(dec) - ta(i) - 12) / 25);
        nb = k >= 150;
        v(nb) = v(nb) + 0.25 * (k(nb) - 150);            % next beat, steeper
        v(nb) = min(v(nb), 0.9);                         % stays below the first peak
        X(r, c, :) = v;
    end
    lat = compute_lat_50(X);
    expect = ta + 5;                                     % 50% of 0 -> 1
    err = abs(lat - expect);
    check('A. no pixel latches the next beat', all(lat(:) < 100), ...
          sprintf('max LAT %.2f frames', max(lat(:))));
    check('A. LAT is the first upstroke crossing', all(err(:) < 1e-9), ...
          sprintf('max error %.2g frames', max(err(:))));

    % ---- B. max dV/dt below the 50% level -----------------------------------
    % 0 -> 0.4 in 2 frames (0.2/frame), then 0.4 -> 1 in 10 frames
    % (0.06/frame), 2-frame plateau at 1, full relaxation afterwards.  The
    % steepest sampled segment always lies inside the foot, and its tangent
    % reaches 0.5 at ta + 2.5 exactly.  (The literal crossing would be
    % ta + 2 + 0.1/0.06.)
    X = zeros([size(ta) nt]);
    for i = 1:numel(ta)
        [r, c] = ind2sub(size(ta), i);
        t = k - ta(i);
        v = zeros(1, nt);
        v(t > 0 & t <= 2)  = 0.2 * t(t > 0 & t <= 2);
        v(t > 2 & t <= 12) = 0.4 + 0.06 * (t(t > 2 & t <= 12) - 2);
        v(t > 12 & t < 14) = 1;
        v(t >= 14) = exp(-(t(t >= 14) - 14) / 15);
        X(r, c, :) = v;
    end
    lat = compute_lat_50(X);
    expect = ta + 2.5;
    err = abs(lat - expect);
    check('B. steep foot extrapolated to 50%', all(err(:) < 1e-9), ...
          sprintf('max error %.2g frames', max(err(:))));

    % ---- D. one-frame seam after the AP (the rabbit 4-18 case) ---------------
    % Every pixel: upstroke 0 -> 1 over 20 frames (0.05/frame), decay to a 0.13
    % floor, a dip to 0.05 on the last frame (so the population crosses 10%
    % only there and the search runs to the end of the window).  One column
    % in four also carries a +0.12 one-frame step at frame 154, steeper than
    % its own upstroke.  Before the fix those pixels got LAT ~156.
    X = zeros([size(ta) nt]);
    seam = false(size(ta));  seam(:, 1:4:end) = true;
    for i = 1:numel(ta)
        [r, c] = ind2sub(size(ta), i);
        v = min(max((k - ta(i)) / 20, 0), 1);            % plateau at 1 for 2 frames
        dec = k >= ta(i) + 22;
        v(dec) = 0.13 + 0.87 * exp(-(k(dec) - ta(i) - 22) / 20);
        v(end) = 0.05;
        if seam(i), v(k >= 154) = v(k >= 154) + 0.12; end
        X(r, c, :) = v;
    end
    lat = compute_lat_50(X);
    expect = ta + 10;
    err = abs(lat - expect);
    check('D. seam pixels keep their own upstroke', all(err(seam) < 1e-9), ...
          sprintf('max error %.2g frames, max LAT %.2f', max(err(seam)), max(lat(seam))));
    check('D. other pixels unaffected', all(err(~seam) < 1e-9), ...
          sprintf('max error %.2g frames', max(err(~seam))));

    % ---- C. ordinary sigmoid upstroke ---------------------------------------
    X = zeros([size(ta) nt]);
    for i = 1:numel(ta)
        [r, c] = ind2sub(size(ta), i);
        t = k - ta(i);
        up = 1 ./ (1 + exp(-t / 1.5));                   % 50% exactly at ta
        v  = up .* exp(-max(t - 15, 0) / 20);            % relaxes to ~0
        X(r, c, :) = v;
    end
    lat = compute_lat_50(X);
    err = abs(lat - ta);
    check('C. sigmoid upstroke stays sub-frame accurate', all(err(:) < 0.05), ...
          sprintf('max error %.3f frames', max(err(:))));

    fprintf('All checks passed.\n');
end


function check(name, cond, detail)
    if cond
        fprintf('   PASS  %-45s (%s)\n', name, detail);
    else
        error('mv_lat_test:failed', 'FAIL  %s (%s)', name, detail);
    end
end
