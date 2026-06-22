function ps_dynamics = update_ps_database(ps_dynamics, ps)
% PSDynamics is a handle class — this is a reference, not a copy.
%
% Records EVERY tracked PS (no validity gating here).  The validity threshold
% is applied afterwards by finalize_ps_validity, so the full lifespan
% distribution is available first — required for the data-driven FDR threshold
% mode, and harmless for the fixed-cycle mode.

% Cache repeated field accesses
ps_bd1      = ps.birthday(1,1);
ps_bd2      = ps.birthday(1,2);
ps_lifespan = ps.lifespan(1,2);

% ps_info: one row per tracked PS (duration, distance, displacement).
% end+1 on an empty matrix correctly appends as row 1 — no isempty guard needed
ps_dynamics.ps_info(end+1,:) = [ps_bd1, ps_bd2, ps_lifespan, ps.distance, ps.displacement];

% all_tracks: keep the full path/wf_path for every PS so validity can be
% decided later against any threshold.
n = size(ps_dynamics.all_tracks, 1);
ps_dynamics.all_tracks{n+1, 1} = [ps_bd1, ps_bd2, ps_lifespan];
ps_dynamics.all_tracks{n+1, 2} = ps.path;
ps_dynamics.all_tracks{n+1, 3} = ps.wf_path;
