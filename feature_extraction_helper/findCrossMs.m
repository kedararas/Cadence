function tMs = findCrossMs(post, level, dt)
%FINDCROSSMS  First time (ms, relative to peak at index 1) at which
% the post-peak signal falls AT OR BELOW `level`. Linear interpolation
% between bracketing samples.
    idx = find(post <= level, 1, 'first');
    if isempty(idx) || idx == 1
        tMs = NaN;
        return;
    end
    % linear interp between idx-1 and idx
    y1 = post(idx-1); y2 = post(idx);
    if y1 == y2
        frac = 0;
    else
        frac = (y1 - level) / (y1 - y2);
    end
    tMs = ((idx - 2) + frac) * dt * 1000;
end