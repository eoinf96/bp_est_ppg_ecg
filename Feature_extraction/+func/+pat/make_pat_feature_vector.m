function feat_tab = make_pat_feature_vector(PAT,t_window_start, t_window_end, configs)
%  This function returns a feature vector of averaged PAT values
%
% INPUT: PPG: This struct must contain pw_inds
%        t_window_start: vector of window start times within which pw inds are averaged.
%        t_window_end: vector of window end times within which pw inds are averaged.
%        configs: configs struct, see below for details
%
% OUTPUT: feat_tab: PAT  feature vector table
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
default_configs.sqi_threshold = 0;
default_configs.window_average_sqi_thresh = 0.4;
configs = func.aux_functions.update_with_default_opts(configs, default_configs);
%% Run loop
num_windows = length(t_window_start);
if ~isfield(PAT, 'good_beat'); PAT.good_beat = PAT.ts_beat; PAT.good_beat(PAT.sqi_beat <=configs.sqi_threshold) = nan;     end

feat_store = cell(num_windows, 1);
for wind_idx = 1:num_windows
    ind_window = and(PAT.t_beat > t_window_start(wind_idx), PAT.t_beat < t_window_end(wind_idx));
    sqi_w = mean(PAT.sqi_beat(ind_window), 'omitnan');
    if  sqi_w >= configs.window_average_sqi_thresh
        feat_store{wind_idx} = mean(PAT.good_beat(ind_window), 1, 'omitnan');
    else
        feat_store{wind_idx} = nan;
    end
    
end

feat_tab = array2table(cell2mat(feat_store), 'VariableNames', {'PAT'});
end

