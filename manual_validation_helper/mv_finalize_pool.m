function manifest = mv_finalize_pool(manifest_file, varargin)
%MV_FINALIZE_POOL  Turn the lead reviewer's accepted candidates into the pool.
%
%   manifest = mv_finalize_pool(manifest_file, 'Marks', lead_marks_file)
%   manifest = mv_finalize_pool(manifest_file, 'Active', logical_vector)
%
%   Called by mv_mark at the end of a 'Lead' session; also callable by hand,
%   e.g. when the lead stopped early and the study will proceed with fewer
%   pixels than the target.
%
%   What it does
%     - accepted  = candidates the lead completed, in manifest order, capped at
%                   manifest.n_target so the pre-registered stopping rule holds
%                   even if the lead marked past it;
%     - active    = those pixels; everything else is out of the study;
%     - grid mode : SNR strata are assigned now, by rank over the ACTIVE pixels
%                   (equal-count groups, so a stratum index means the same thing
%                   as in stratified mode);
%     - records the lead's identity, marks file, skipped candidates and the
%                   number visited, then sets pool_status = 'final' and saves.
%
%   A final pool is not re-opened: re-running on it errors unless 'Force' is
%   set, because other reviewers may already be marking it.
%
%   Name-value
%     'Marks'   the lead's marks file (normal use)
%     'Active'  explicit logical vector over candidates (used by
%               mv_sample_pixels when Candidates = 1, and for tests)
%     'Force'   re-finalise an already-final pool (default false)
%     'Quiet'   suppress the summary (default false)

    p = inputParser;
    p.addParameter('Marks',  '',    @(v) ischar(v) || isstring(v));
    p.addParameter('Active', [],    @(v) isempty(v) || islogical(v));
    p.addParameter('Force',  false, @islogical);
    p.addParameter('Quiet',  false, @islogical);
    p.parse(varargin{:});
    o = p.Results;

    manifest_file = char(manifest_file);
    M = load(manifest_file);
    if ~isfield(M, 'manifest')
        error('mv_finalize_pool:badManifest', '%s does not contain a ''manifest'' struct.', manifest_file);
    end
    manifest = M.manifest;
    n = size(manifest.pixels, 1);

    if ~isfield(manifest, 'pool_status')
        error('mv_finalize_pool:oldManifest', ...
              '%s predates candidate pools (no pool_status). Re-run mv_sample_pixels.', manifest_file);
    end
    if strcmp(manifest.pool_status, 'final') && ~o.Force
        error('mv_finalize_pool:alreadyFinal', ...
              ['%s is already final (%d pixels, finalised %s by %s). Other reviewers may be ' ...
               'marking it. Pass ''Force'', true only if you mean to redefine the pool.'], ...
              manifest_file, sum(manifest.active), manifest.pool_finalized, manifest.lead_reviewer);
    end

    if ~isempty(o.Active)
        if numel(o.Active) ~= n
            error('mv_finalize_pool:activeSize', ...
                  '''Active'' has %d entries; the manifest has %d candidates.', numel(o.Active), n);
        end
        active  = o.Active(:);
        skipped = zeros(0, 1);
        visited = n;
        lead    = '';
        lead_file = '';
    elseif ~isempty(o.Marks)
        lead_file = char(o.Marks);
        R = load(lead_file);
        if ~isfield(R, 'marks')
            error('mv_finalize_pool:badMarks', '%s contains no ''marks'' struct.', lead_file);
        end
        m = R.marks;
        if ~same_file(m.manifest_file, manifest_file)
            error('mv_finalize_pool:manifestMismatch', ...
                  'Marks file %s references manifest\n  %s\nbut you passed\n  %s', ...
                  lead_file, m.manifest_file, manifest_file);
        end
        if numel(m.done) ~= n
            error('mv_finalize_pool:sizeMismatch', ...
                  'Marks file covers %d pixels; manifest has %d candidates.', numel(m.done), n);
        end
        accepted = find(m.done(:));                  % manifest order
        if numel(accepted) > manifest.n_target
            if ~o.Quiet
                fprintf('Lead accepted %d, target is %d: keeping the first %d in manifest order.\n', ...
                        numel(accepted), manifest.n_target, manifest.n_target);
            end
            accepted = accepted(1:manifest.n_target);
        end
        active = false(n, 1);
        active(accepted) = true;
        skipped = find(m.skipped(:));
        visited = max([accepted; skipped; 0]);
        lead    = char(m.reviewer);
    else
        error('mv_finalize_pool:noInput', 'Pass either ''Marks'' or ''Active''.');
    end

    n_active = sum(active);
    if n_active == 0
        error('mv_finalize_pool:empty', 'No accepted pixels — nothing to finalise.');
    end
    if n_active < manifest.n_target
        warning('mv_finalize_pool:short', ...
                ['Pool has %d pixels; the target was %d. The candidate list is exhausted ' ...
                 '(or the lead stopped early). Proceed with %d, or re-sample with a larger ' ...
                 '''Candidates'' factor BEFORE any other reviewer starts.'], ...
                n_active, manifest.n_target, n_active);
    end

    % Post-stratify the grid pool by rank.  See mv_sample_pixels for why rank
    % rather than percentile value (ties would empty a stratum).
    if strcmp(manifest.sampling, 'grid')
        K = manifest.num_strata;
        if n_active < K
            error('mv_finalize_pool:tooFewForStrata', ...
                  'Pool has %d pixels, fewer than the %d requested strata.', n_active, K);
        end
        idx        = find(active);
        s_sel      = manifest.snr(idx);
        [~, ord]   = sort(s_sel, 'ascend');
        cuts       = round(linspace(0, n_active, K + 1));
        stratum    = nan(n, 1);
        for k = 1:K
            stratum(idx(ord(cuts(k)+1 : cuts(k+1)))) = k;
        end
        manifest.stratum   = stratum;
        manifest.snr_edges = [-inf, arrayfun(@(k) s_sel(ord(cuts(k+1))), 1:K-1), inf];
    end

    manifest.active          = active;
    manifest.pool_status     = 'final';
    manifest.lead_reviewer   = lead;
    manifest.lead_marks_file = lead_file;
    manifest.lead_skipped    = skipped;
    manifest.lead_visited    = visited;
    manifest.pool_finalized  = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

    save(manifest_file, 'manifest');
    manifest.manifest_file = manifest_file;

    if ~o.Quiet
        fprintf('Pool finalised: %d of %d target', n_active, manifest.n_target);
        if ~isempty(lead)
            fprintf(' (lead %s: %d candidates visited, %d skipped = %.0f%% unmarkable)', ...
                    lead, visited, numel(skipped), 100 * numel(skipped) / max(visited, 1));
        end
        fprintf('\n');
        if strcmp(manifest.sampling, 'grid')
            cnt = arrayfun(@(t) sum(active & manifest.tier == t), 1:max(manifest.tier));
            fprintf('  tiers in pool : %s\n', mat2str(cnt));
        end
        fprintf('  SNR in pool   : %.1f-%.1f, %d strata\n', ...
                min(manifest.snr(active)), max(manifest.snr(active)), manifest.num_strata);
        fprintf('  manifest      : %s\n', manifest_file);
        fprintf('Other reviewers can now run mv_mark on this manifest.\n');
    end
end


function tf = same_file(a, b)
    a = char(a); b = char(b);
    tf = strcmp(a, b);
    if ~tf
        [~, na, ea] = fileparts(a); [~, nb, eb] = fileparts(b);
        tf = strcmp([na ea], [nb eb]);
    end
end
