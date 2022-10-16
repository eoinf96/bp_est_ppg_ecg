function feat_tab = make_ECG_feature_vector(ECG,t_window_start, t_window_end, opts)
if nargin <4
    opts = struct();
end
default_opts.scale_list = 2:2:8;
default_opts.do_HRV_feats = 0;
default_opts.HRV_list  = {'SDRR', 'pRR50', 'RMSSD', 'poincare_SDD1', 'poincare_SDD2', 'poincare_area', 'poincare_SDD1_SDD2'}; % only relevant if opts.do_HRV_feats ==1
opts = func.aux_functions.update_with_default_opts(opts, default_opts);
%% Run loop
num_windows = length(t_window_start);

feat_store = cell(num_windows, 1);
for wind_idx = 1:num_windows
    ind_t = and(ECG.t > t_window_start(wind_idx), ECG.t < t_window_end(wind_idx));
    [feat_store{wind_idx}, names] = func.ecg.get_ECG_features(ECG.ts(ind_t), opts);
end

feat_tab = array2table(cell2mat(feat_store')', 'VariableNames', names);
if opts.do_HRV_feats
    RR.ts = ECG.t(ECG.peaks(2:end)) - ECG.t(ECG.peaks(1:end-1));
    RR.t = ECG.t(ECG.peaks(2:end));

    HRV_feats = func.HRV.get_HRV_features(RR.ts, RR.t, t_window_start, t_window_end);
    HRV_feats = struct2table(HRV_feats);
    HRV_feats = HRV_feats(:, opts.HRV_list);
    feat_tab= [feat_tab ,HRV_feats];
end



end