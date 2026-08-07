classdef WavefrontDynamics < handle
    % WavefrontDynamics  Handle class accumulating wavefront tracking results.
    %
    % Using handle semantics means update_wf_database modifies this object
    % in-place — no struct copy on every call as the database grows.

    properties
        wf_count
        df
        frame_rate
        pixel_size        = []   % mm/px, if known; enables physical reentry tolerances
        wf_size_duration  = []
        wf_path           = {}
        wf_fractionations = []
        wf_collisions     = []
        wf_blocks         = []
        wf_breakthroughs  = []
        wf_multiplicity   = []
        wf_repeatability  = []
        wf_reentry        = []
    end

    methods
        function obj = WavefrontDynamics(wf_count, df, frame_rate, pixel_size)
            obj.wf_count    = wf_count;
            obj.df          = df;
            obj.frame_rate  = frame_rate;
            if nargin >= 4
                obj.pixel_size = pixel_size;
            end
        end
    end
end
