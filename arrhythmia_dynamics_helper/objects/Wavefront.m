classdef Wavefront < handle
    % Wavefront  Handle class representing a single tracked wavefront.
    %
    % Using handle semantics means all references to the same Wavefront
    % object share state — no copy-on-write overhead when passing between
    % functions or storing in cell arrays. Modifications to a retrieved
    % reference are immediately visible everywhere that reference is held.

    properties
        parent
        child
        mate
        fling
        neighbors
        lifespan
        frame
        index
        length
        location
        path
        birthday
        cause_of_death
    end

    methods
        function obj = Wavefront(frame, index)
            obj.frame = frame;
            obj.index = index;
        end
    end
end

