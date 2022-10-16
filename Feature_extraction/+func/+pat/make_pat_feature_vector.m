function feat_tab = make_pat_feature_vector(PAT,t_window_start, t_window_end, opts)
    narginchk(3, inf)
    if nargin < 4
       opts = struct();
    end
    default_opts.sqi_threshold = 0;
    opts = func.aux_functions.update_with_default_opts(opts, default_opts);
    %% Run loop
    num_windows = length(t_window_start);
    if ~isfield(PAT, 'good_beat'); PAT.good_beat = PAT.ts_beat; PAT.good_beat(PAT.sqi_beat <=opts.sqi_threshold) = nan;     end

    feat_store = cell(num_windows, 1);
    for wind_idx = 1:num_windows
        ind_window = and(PAT.t_beat > t_window_start(wind_idx), PAT.t_beat < t_window_end(wind_idx));
        feat_store{wind_idx} = mean(PAT.good_beat(ind_window), 1, 'omitnan');
    end

    feat_tab = array2table(cell2mat(feat_store), 'VariableNames', {'PAT'});
end

