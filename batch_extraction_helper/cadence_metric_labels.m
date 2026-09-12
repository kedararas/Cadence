function L = cadence_metric_labels()
%CADENCE_METRIC_LABELS  Master list of per-recording summary metrics.
%
%   L = cadence_metric_labels()  ->  N x 3 cell {name, label, description}
%
%   name   : valid MATLAB identifier used as a table variable / struct field
%   label  : spreadsheet column header (unit in square brackets)
%   desc   : how the number is obtained from the *-metrics ep_metrics struct
%
%   Every value is the MEDIAN over finite pixels of the masked (slot 1)
%   map, unless the description says "scalar".  Labels that were already
%   used by the previous Python summary (extract_medians.py) are kept
%   verbatim so downstream scripts keep working:
%     Activation time (V) [ms], APD (V) [ms], Repolarization time (V) [ms],
%     AP rise time (V) [ms], Voltage-Ca delay [ms], Ca transient duration [ms],
%     Ca decay time [ms], Ca rise time [ms], Ca decay tau [ms],
%     Dominant frequency (V) [Hz], Regularity index (V), Organization index (V),
%     APD spatial gradient (V), APD restitution slope (V),
%     APD alternans ratio (V), AP amplitude ratio (V), Alternans phase (V).

    V  = 'voltage camera (CAM1)';
    Ca = 'calcium camera (CAM2)';
    L = { ...
    % ---- voltage ------------------------------------------------------------
    'act_time_v',        'Activation time (V) [ms]',             ['act_times: LAT at 50% upstroke of the ensemble beat, ' V];
    'rep_time_v',        'Repolarization time (V) [ms]',         'rep_data: time of 80% repolarization from window start';
    'ap_rise_v',         'AP rise time (V) [ms]',                'ap_rise_times: 20% -> 90% upstroke time';
    'apd_v',             'APD (V) [ms]',                         'apd_data: APD80 (repolarization minus activation)';
    'apd50_v',           'APD50 (V) [ms]',                       'apd50_data: APD at 50% repolarization, same engine as APD80 (computed on the fly from CAM1_average if the field is absent)';
    'cv_v',              'Conduction velocity (V) [cm/s]',       'local_cv{1,1}: inverse-gradient local CV map saved by Feature Extraction, interior only (8 px erosion, as in Signal Analysis)';
    'cv_polyfit_v',      'Conduction velocity polyfit (V) [cm/s]', 'conduction_velocity(act map, ''polyfit'') recomputed from the saved activation map, interior only (the Conduction Velocity module default engine)';
    % ---- calcium ------------------------------------------------------------
    'vc_delay',          'Voltage-Ca delay [ms]',                ['vc_delay: Ca activation minus Vm activation (signed), ' Ca];
    'ca_decay_ca',       'Ca decay time [ms]',                   'ca_rep_data: time of 80% decay from window start';
    'catd_ca',           'Ca transient duration [ms]',           'ca_data: CaTD80';
    'catd50_ca',         'CaTD50 (Ca) [ms]',                     'ca50_data: Ca transient duration at 50% decay (computed on the fly from CAM2_average if the field is absent)';
    'ca_tau_ca',         'Ca decay tau [ms]',                    'ca_tau{1}: single-exponential decay constant fitted 90% -> 20%';
    'ca_tau_1e_ca',      'Ca decay tau 1/e [ms]',                'ca_tau{6}: model-free time from peak to 1/e of amplitude';
    'ca_rise_ca',        'Ca rise time [ms]',                    'ca_rise_times: 20% -> 90% upstroke time';
    % ---- voltage spectral --------------------------------------------------
    'df_v',              'Dominant frequency (V) [Hz]',          'complexity_data{1,1}{1}: pixel-wise dominant frequency';
    'ri_v',              'Regularity index (V)',                 'complexity_data{1,1}{2}: regularity index (0..1)';
    'oi_v',              'Organization index (V)',               'complexity_data{1,1}{3}: organization index (0..1)';
    % ---- voltage alternans -------------------------------------------------
    'apd_alt_ratio_v',   'APD alternans ratio (V)',              'substrate.alt_ratio_map: |APD80 odd-even| / mean APD80';
    'apd80_alt_ratio_v', 'APD80 alternans ratio map (V)',        'analyzeAPAlternans APD80_ratio_map';
    'apd50_alt_ratio_v', 'APD50 alternans ratio (V)',            'analyzeAPAlternans APD50_ratio_map: |APD50 odd-even| / mean APD50';
    'amp_alt_ratio_v',   'AP amplitude ratio (V)',               'analyzeAPAlternans amp_ratio_map: |amplitude odd-even| / mean amplitude';
    'dvdt_alt_ratio_v',  'dV/dt alternans ratio (V)',            'analyzeAPAlternans dvdt_ratio_map';
    'sai_v',             'Spectral alternans index (V)',         'analyzeAPAlternans spectral_map: power at 0.5 cycles/beat of the per-beat APD series';
    'alt_phase_v',       'Alternans phase (V)',                  'substrate.phase_map: +1 odd longer, -1 even longer (median over pixels)';
    'alt_sig_frac_v',    'Alternans significant fraction (V)',   'scalar: fraction of tissue pixels with APD80 ratio > 0.10 and paired-t p < 0.05';
    'concordance_v',     'Alternans concordance ratio (V)',      'scalar substrate.concordance_ratio: fraction of significant pixels in the majority phase (< 0.80 = discordant)';
    'apd_gradient_v',    'APD spatial gradient (V)',             'substrate.gradient_mag: Sobel gradient of the mean APD map (ms/px)';
    'restitution_v',     'APD restitution slope (V)',            'substrate.restitution_map: local dAPD/dDI from odd/even beat pairs';
    % ---- calcium alternans -------------------------------------------------
    'ca_amp_alt_ratio',  'CaT amplitude alternans ratio (Ca)',   'analyzeCaTransientAlternans amp_ratio_map';
    'ca_d50_alt_ratio',  'CaTD50 alternans ratio (Ca)',          'analyzeCaTransientAlternans D50_ratio_map';
    'ca_d80_alt_ratio',  'CaTD80 alternans ratio (Ca)',          'analyzeCaTransientAlternans D80_ratio_map (check level coverage at fast pacing)';
    'ca_release_alt',    'Ca release alternans (Ca)',            'analyzeCaTransientAlternans release_alt_map: 1 - S/L (Wang 2014)';
    'ca_load_alt',       'Ca load alternans (Ca)',               'analyzeCaTransientAlternans load_alt_map: D/L (Wang 2014)';
    'ca_sai',            'Spectral alternans index (Ca)',        'analyzeCaTransientAlternans spectral_map';
    'ca_alt_sig_frac',   'Alternans significant fraction (Ca)',  'scalar: fraction of tissue pixels with CaT amplitude ratio > 0.10 and paired-t p < 0.05';
    'ca_alt_phase',      'Alternans phase (Ca)',                 'analyzeCaTransientAlternans phase_map (median over pixels)';
    'ca_ap_coupling',    'Ca-AP coupling r',                     'substrate.coupling_map: per-pixel Pearson r between per-beat Ca amplitude and APD';
    'ca_ap_inphase',     'Ca-AP in-phase fraction',              'scalar substrate.coupling_fraction';
    % ---- composite ---------------------------------------------------------
    'risk_map',          'Arrhythmia risk (composite)',          'substrate.risk_map (combined Vm + Ca substrate), median over pixels';
    'risk_global',       'Arrhythmia risk global',               'scalar substrate.risk_global';
    % ---- beat-windowed (per-beat maps from the alternans analysis, median over beats then pixels) ----
    'apd50_bw_v',        'APD50 (V, beat-windowed) [ms]',        'analyzeAPAlternans APD_beat at 50%: each beat in its own stimulus-to-stimulus window, LAT 50% upstroke to 50% repolarization; per-pixel median over beats';
    'apd80_bw_v',        'APD80 (V, beat-windowed) [ms]',        'analyzeAPAlternans APD_beat at 80%';
    'ca_ttp_bw',         'Ca time-to-peak (beat-windowed) [ms]', 'analyzeCaTransientAlternans TTP_beat: beat trough (onset) to peak';
    'ca_d50_bw',         'Ca decay 50% from peak (beat-windowed) [ms]', 'analyzeCaTransientAlternans CaTD_beat at 50%: peak to 50% decay toward the beat''s own trough';
    'ca_d80_bw',         'Ca decay 80% from peak (beat-windowed) [ms]', 'analyzeCaTransientAlternans CaTD_beat at 80%';
    'catd50_bw',         'CaTD50 (beat-windowed) [ms]',          'TTP + decay-50 per beat and pixel: onset to 50% decay, comparable to the ensemble CaTD50';
    'catd80_bw',         'CaTD80 (beat-windowed) [ms]',          'TTP + decay-80 per beat and pixel: onset to 80% decay, comparable to the ensemble CaTD80';
    'ca_diast_bw',       'Ca diastolic level (beat-windowed)',   'analyzeCaTransientAlternans diastolic_beat: per-beat trough level in the conditioned signal units; rises when the transient does not relax within the cycle';
    'ca_ensemble_valid', 'Ca ensemble metrics valid (flag)',     'scalar: 1 when the ensemble-window Ca metrics pass every test: the app guard on the averaged beat (check_analysis_window: not wrapped, decay fraction >= 0.5), V-Ca delay >= 0, CaTD valid fraction >= 0.5, CaTD80 >= 15 ms, and within 50% of the beat-windowed CaTD80 when available; 0 = treat Ca decay time / CaTD / Ca rise / tau / V-Ca delay as invalid';
    'v_ensemble_valid',  'V ensemble metrics valid (flag)',      'scalar: 1 when the ensemble-window V metrics pass every test: the app guard on the averaged beat, APD valid fraction >= 0.5, APD80 >= 10 ms, and within 50% of the beat-windowed APD80 when available';
    % ---- QC ----------------------------------------------------------------
    'n_beats',           'Beats analysed (alternans)',           'scalar: number of beats segmented from analog1 for the alternans analysis';
    'tissue_frac_v',     'Tissue fraction (V mask)',             'scalar: fraction of the frame inside the adaptive SNR mask, voltage camera';
    'tissue_frac_ca',    'Tissue fraction (Ca mask)',            'scalar: fraction of the frame inside the adaptive SNR mask, calcium camera';
    'apd_valid_frac_v',  'APD valid pixel fraction (V)',         'scalar: pixels with a finite APD80 / pixels with a finite LAT';
    'catd_valid_frac_ca','CaTD valid pixel fraction (Ca)',       'scalar: pixels with a finite CaTD80 / pixels inside the Ca mask';
    'v_polarity_conf',   'V polarity confidence (stimulus)',    'scalar: post-stimulus deflection of the voltage camera (check_polarity_paced), +1 upright to -1 inverted. NEGATIVE means the stored conditioned data is upside-down and every voltage metric for that recording is invalid; re-condition it';
    'ca_polarity_conf',  'Ca polarity confidence (stimulus)',   'scalar: same test on the calcium camera. Conditioning decides polarity per recording from signal shape, which fails at fast pacing because the transient fills the cycle; a negative value here marks a recording it flipped in error';
    'pacing_hz',         'Pacing rate measured [Hz]',            'scalar: 1000 / median inter-stimulus interval of analog1 (ms)';
    'capture_ratio',     'Capture ratio (DF / pacing rate)',     'scalar: voltage dominant frequency / measured pacing rate; ~1 = 1:1 capture, ~0.5 = 2:1 block. Functional refractory period = shortest CL still at ~1';
    'v_ca_identical',    'V and Ca maps identical (QC flag)',    'scalar: 1 if the APD80 map equals the CaTD80 map or the V and Ca activation maps are equal (would mean both cameras were read from the same data), else 0';
    };
end
