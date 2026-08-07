function recs = validate_metrics(d, recording, species)
%VALIDATE_METRICS  QC gate for the feature-extraction stage output.
%
%   recs = validate_metrics(cmos_all_data, recording)
%   recs = validate_metrics(cmos_all_data, recording, species)
%
%   Run this in the batch loop right after feature/arrhythmia extraction has
%   populated cmos_all_data (df_map, ps_dynamics, pixel_size, ...), before it
%   is saved.  Returns an array of qc_check records.
%
%   species (optional) one of "rat","rabbit","pig","human" — only used to
%   widen/narrow the expected dominant-frequency band in the message; the
%   pass/fail band stays physiological and generous so it never hard-fails a
%   real but unusual recording.
%
%   Checks
%     pixel_size     present and physically plausible (mm/px).
%     df_rotor       p75 of the DF map (the rotor-frequency estimate) finite
%                    and inside a physiological band.
%     rotor_tracking PS were tracked at all.
%     rotor_valid    distinguishes "no sustained rotors found" (WARN, a real
%                    negative result) from "tracking produced nothing" (FAIL).
%
%   See qc_check for the record shape and status semantics.

    if nargin < 3, species = ""; end
    species = lower(string(species));

    recs = repmat(qc_check("","","" ,"PASS",NaN,""), 0, 1);
    name = string(recording);

    function add(check, status, value, message)
        recs(end+1, 1) = qc_check(name, "features", check, status, value, message);
    end

    % ---- pixel size (mm/px) ----
    if isfield(d, 'pixel_size') && isscalar(d.pixel_size) && isfinite(d.pixel_size) && d.pixel_size > 0
        px = double(d.pixel_size);
        if px > 0.005 && px < 1.0           % 5 um .. 1 mm per pixel
            add("pixel_size", "PASS", px, sprintf("%.4f mm/px", px));
        else
            add("pixel_size", "WARN", px, sprintf("%.4f mm/px implausible (check FOV)", px));
        end
    else
        add("pixel_size", "FAIL", NaN, "pixel_size missing/non-positive (thresholds fall back to geometry)");
    end

    % ---- rotor frequency from DF map ----
    if isfield(d, 'df_map') && ~isempty(d.df_map) && isnumeric(d.df_map)
        df_rotor = prctile(d.df_map(:), 75);   % NaNs ignored
        % Expected DF bands (approx) by species, for the message only.
        band = struct('rat',[8 35], 'rabbit',[5 25], 'pig',[3 15], 'human',[3 12]);
        exptxt = "";
        if species ~= "" && isfield(band, species)
            b = band.(species);
            exptxt = sprintf(" (expected ~%d-%d Hz for %s)", b(1), b(2), species);
        end
        if ~isfinite(df_rotor) || df_rotor <= 0
            add("df_rotor", "FAIL", df_rotor, "DF map empty / non-positive");
        elseif df_rotor < 1 || df_rotor > 50
            add("df_rotor", "WARN", df_rotor, ...
                sprintf("%.1f Hz outside physiological [1,50] Hz%s", df_rotor, exptxt));
        else
            add("df_rotor", "PASS", df_rotor, sprintf("%.1f Hz%s", df_rotor, exptxt));
        end
    else
        add("df_rotor", "FAIL", NaN, "df_map missing");
    end

    % ---- rotor tracking / validity ----
    if isfield(d, 'ps_dynamics') && ~isempty(d.ps_dynamics)
        psd      = d.ps_dynamics;
        n_tracks = size(psd.ps_info, 1);
        n_valid  = size(psd.path, 1);
        if n_tracks == 0
            add("rotor_tracking", "FAIL", 0, "no phase singularities tracked");
        elseif n_valid == 0
            % Real negative result, not a failure: PS existed but none met the
            % validity threshold (no sustained rotors).
            add("rotor_valid", "WARN", 0, ...
                sprintf("%d PS tracked, 0 met validity threshold (no sustained rotors)", n_tracks));
        else
            add("rotor_valid", "PASS", n_valid, ...
                sprintf("%d valid rotor(s) of %d tracked PS", n_valid, n_tracks));
        end
    else
        add("rotor_tracking", "FAIL", NaN, "ps_dynamics missing");
    end

    if isempty(recs)
        add("structure", "FAIL", NaN, "no recognised metric fields");
    end
end
