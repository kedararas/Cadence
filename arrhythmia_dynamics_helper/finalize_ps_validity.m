function ps_dynamics = finalize_ps_validity(ps_dynamics, valid_ps_threshold)
%FINALIZE_PS_VALIDITY  Populate path/valid_ps from all_tracks given a threshold.
%
%   ps_dynamics = finalize_ps_validity(ps_dynamics, valid_ps_threshold)
%
%   Runs after tracking.  Every tracked PS is stored in ps_dynamics.all_tracks;
%   this pass keeps only those whose lifespan >= valid_ps_threshold, writing
%   them into ps_dynamics.path and tagging ps_dynamics.valid_ps along each
%   surviving track.  Idempotent — path and valid_ps are reset first, so it can
%   be re-run with a different threshold without re-tracking.

    % Reset outputs so repeated calls are clean.
    ps_dynamics.path = {};
    for f = 1:numel(ps_dynamics.valid_ps)
        if ~isempty(ps_dynamics.valid_ps{f})
            ps_dynamics.valid_ps{f}(:) = 0;
        end
    end

    for t = 1:size(ps_dynamics.all_tracks, 1)
        meta     = ps_dynamics.all_tracks{t, 1};   % [bd_frame, bd_index, lifespan]
        pth      = ps_dynamics.all_tracks{t, 2};
        wf_path  = ps_dynamics.all_tracks{t, 3};
        lifespan = meta(3);

        if lifespan >= valid_ps_threshold
            n = size(ps_dynamics.path, 1);
            ps_dynamics.path{n+1, 1} = meta;
            ps_dynamics.path{n+1, 2} = pth;
            ps_dynamics.path{n+1, 3} = wf_path;

            for i = 1:size(pth, 1)
                frame = pth(i, 1);
                index = pth(i, 2);
                ps_dynamics.valid_ps{frame}(index, 1) = 1;
            end
        end
    end
end
