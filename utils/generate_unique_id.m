function unique_id = generate_unique_id
%GENERATE_UNIQUE_ID  Create a globally-unique identifier string.
%
%   unique_id = generate_unique_id
%
%   Returns a random RFC-4122 UUID (via Java's java.util.UUID) as a
%   character row vector, e.g. '3f2504e0-4f89-41d3-9a0c-0305e82c3301'.
%   Used to tag interactively-created UI objects in the CADENCE modules
%   (signal-window ROIs, user data masks, CV point selections) so each can
%   be referenced, resized, or deleted independently. Collisions are
%   astronomically unlikely, so no registry of issued IDs is needed.
%
%   Inputs
%     (none)
%
%   Output
%     unique_id  1xN char, the UUID in canonical 8-4-4-4-12 hex form.

    temp = java.util.UUID.randomUUID;
    unique_id = char(temp.toString);

end