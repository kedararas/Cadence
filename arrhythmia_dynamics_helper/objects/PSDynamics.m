classdef PSDynamics < handle
    % PSDynamics  Handle class accumulating phase singularity tracking results.
    %
    % Handle semantics mean update_ps_database modifies this object in-place —
    % no struct copy on every call as the database grows.

    properties
        ps_info    = []
        path       = {}
        valid_ps   = {}
        all_tracks = {}   % every tracked PS regardless of validity:
                          %   {n,1}=[birthday_frame, birthday_index, lifespan]
                          %   {n,2}=path  {n,3}=wf_path
                          % Lets the validity threshold be applied AFTER tracking
                          % (needed for the data-driven FDR threshold mode).
    end

    methods
        function obj = PSDynamics(valid_ps)
            obj.valid_ps = valid_ps;
        end
    end
end
