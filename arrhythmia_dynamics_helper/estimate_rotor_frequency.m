function df_rotor = estimate_rotor_frequency(df_map, df_percentile)
%ESTIMATE_ROTOR_FREQUENCY  Rotor rotation frequency (Hz) from a DF map.
%
%   df_rotor = estimate_rotor_frequency(df_map)
%   df_rotor = estimate_rotor_frequency(df_map, df_percentile)
%
%   Cardiac optical-mapping DF maps are typically BIMODAL: a low-frequency
%   far-field cluster and a faster cluster where the reentrant cores sit.
%   The plain median is dragged into the slow cluster by the far-field
%   pixels, which inflates the rotation period (and any threshold derived
%   from it) several-fold and makes real rotors fail validity checks.
%
%   Taking a high percentile of the DF map places the estimate in the fast,
%   rotor-bearing tissue instead.  This is operator-independent (no pixel
%   click) and batch-safe.
%
%   Inputs
%     df_map         dominant-frequency map (Hz), any size; NaNs ignored.
%     df_percentile  percentile to use (default 75).  Raise toward 90 if
%                    slow far-field pixels still dominate; lower it if
%                    transients are being admitted as valid rotors.
%
%   Output
%     df_rotor       scalar rotor frequency (Hz).

    if nargin < 2 || isempty(df_percentile)
        df_percentile = 75;
    end

    df_rotor = prctile(df_map(:), df_percentile);   % NaNs ignored by prctile

    % Guard against a degenerate/empty map producing NaN or 0 Hz, which would
    % make any frame threshold Inf/NaN and silently reject every rotor.
    if ~isfinite(df_rotor) || df_rotor <= 0
        df_rotor = median(df_map(:), 'omitnan');
    end
end
