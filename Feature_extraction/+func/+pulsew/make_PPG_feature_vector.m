function feat_tab = make_PPG_feature_vector(PPG,t_window_start, t_window_end, configs)
%  This function returns a feature vector formed from PPG indices typically
%  used in BP estimation
%
% INPUT: PPG: This struct must contain pw_inds
%        t_window_start: vector of window start times within which pw inds are averaged.
%        t_window_end: vector of window end times within which pw inds are averaged.
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
narginchk(3, inf)
if nargin < 4
    configs = struct();
end
default_configs.do_norm = 1;
default_configs.window_average_sqi_thresh = 0.4;
configs = func.aux_functions.update_with_default_opts(configs, default_configs);
%% Error check
if configs.do_norm
    if ~isfield(PPG, 'norm_pw_inds')
        PPG.norm_pw_inds = func.pulsew.get_ppg_indices(PPG, 1);
    end
    feats = PPG.norm_pw_inds;
else
    if ~isfield(PPG, 'pw_inds')
        PPG.pw_inds = func.pulsew.get_ppg_indices(PPG, 0);
    end
    feats = PPG.pw_inds;
end

if ~isfield(PPG, 't_beat')
    PPG.t_beat = PPG.t(PPG.peaks);
end

%% Run loop
feat_names = fieldnames(feats);
feats = struct2array(feats);
num_windows = length(t_window_start);

feat_store = cell(num_windows, 1);
for wind_idx = 1:num_windows
    ind_window = and(PPG.t_beat > t_window_start(wind_idx), PPG.t_beat < t_window_end(wind_idx));
    sqi_w = mean(PPG.sqi_beat(ind_window), 'omitnan');
    if  sqi_w >= configs.window_average_sqi_thresh
        feat_store{wind_idx} = mean(feats(ind_window, :), 1, 'omitnan');
    else
        feat_store{wind_idx} = nan(1, length(feat_names));
    end
end

feat_tab = array2table(cell2mat(feat_store), 'VariableNames', feat_names);
end

