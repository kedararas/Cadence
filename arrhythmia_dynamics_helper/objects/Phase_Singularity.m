classdef Phase_Singularity < handle
    % Phase_Singularity  Handle class representing a single tracked phase singularity.
    %
    % Handle semantics eliminate copy-on-write overhead when passing between
    % functions or storing in cell arrays.

    properties
        wavefront
        lifespan
        frame
        index
        location
        displacement
        distance
        path
        wf_path
        charge
        birthday
    end

    methods
        function obj = Phase_Singularity(frame, index)
            obj.frame = frame;
            obj.index = index;
        end
    end
end

