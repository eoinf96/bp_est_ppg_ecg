function feat_tab = make_ECG_feature_vector(ECG,t_window_start, t_window_end, configs)
%  This function returns a feature vector formed from ECG features typically
%  used in BP estimation
%
% INPUT: ECG: ECG signal
%        t_window_start: vector of window start times within which ECG features are computed.
%        t_window_end: vector of window end times within which ECG features are computed.
%        configs: configs struct, see below for details
%
% OUTPUT: feat_tab: Feature vector table
% ---
% Features from the photoplethysmogram and the electrocardiogram for estimating changes in blood pressure.
%
% Released under the GNU General Public License
%
% Copyright (C) 2022  Eoin Finnegan
% University of Oxford, Insitute of Biomedical Engineering, CIBIM Lab
% eoin.finnegan@eng.ox.ac.uk
%
% Referencing this work
%
% Finnegan, E., Davidson, S., Harford, M., Jorge, J., Watkinson, P., Tarassenko, L. and Villarroel, M., 2022. Features from the photoplethysmogram and the electrocardiogram for estimating changes in blood pressure. Submitted to Scientific reports
%
% Relevant literature:
% - Monika Simjanoska et al. “Non-invasive blood pressure estimation from ECG using machine learning techniques”. In: Sensors (Switzerland) 18.4 (2018), p. 1160.
% - Sen Yang et al. “Blood pressure estimation with complexity features from electrocardiogram and photoplethysmogram signals”. In: Optical and Quantum Electronics 52.3 (2020), p. 135.
% - 
narginchk(3, inf)
if nargin <4
    configs = struct();
end
default_configs.MSE_scale_list = 2:2:8; % The different scales to compute mult-scale entropy
default_configs.do_HRV_feats = 1; % Flag whether to return heart rate variability features
default_configs.HRV_list  = {'SDRR', 'pRR50', 'RMSSD', 'poincare_SDD1', 'poincare_SDD2', 'poincare_area', 'poincare_SDD1_SDD2'}; % only relevant if opts.do_HRV_feats ==1
configs = func.aux_functions.update_with_default_opts(configs, default_configs);
%% Run loop
num_windows = length(t_window_start);

feat_store = cell(num_windows, 1);
for wind_idx = 1:num_windows
    ind_t = and(ECG.t > t_window_start(wind_idx), ECG.t < t_window_end(wind_idx));
    [feat_store{wind_idx}, names] = func.ecg.get_ECG_features(ECG.ts(ind_t), configs);
end

feat_tab = array2table(cell2mat(feat_store')', 'VariableNames', names);
if configs.do_HRV_feats
    RR.ts = ECG.t(ECG.peaks(2:end)) - ECG.t(ECG.peaks(1:end-1));
    RR.t = ECG.t(ECG.peaks(2:end));

    HRV_feats = func.HRV.get_HRV_features(RR.ts, RR.t, t_window_start, t_window_end);
    HRV_feats = struct2table(HRV_feats);
    HRV_feats = HRV_feats(:, configs.HRV_list);
    feat_tab= [feat_tab ,HRV_feats];
end



end