function feat_tab = make_PPG_feature_vector(PPG,t_window_start, t_window_end, opts)
    narginchk(3, inf)
    if nargin < 4
       opts = struct();
    end
    default_opts.do_norm = 1;
    opts = func.aux_functions.update_with_default_opts(opts, default_opts);
    %% Run loop
    if opts.do_norm
       feats = PPG.norm_pw_inds; 
    else
       feats = PPG.pw_inds;
    end
    feat_names = fieldnames(feats);
    feats = struct2array(feats);
    num_windows = length(t_window_start);

    feat_store = cell(num_windows, 1);
    for wind_idx = 1:num_windows
        ind_window = and(PPG.t(PPG.peaks) > t_window_start(wind_idx), PPG.t(PPG.peaks) < t_window_end(wind_idx));
        feat_store{wind_idx} = mean(feats(ind_window, :), 1, 'omitnan');
    end

    feat_tab = array2table(cell2mat(feat_store), 'VariableNames', feat_names);
end

