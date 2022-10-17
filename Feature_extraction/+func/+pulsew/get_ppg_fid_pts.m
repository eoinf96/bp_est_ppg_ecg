function [varargout ] = get_ppg_fid_pts(PPG,config, plot_flag)
% This function locates the characteristic fiducial points of the PPG
% Code adapted from P.Charlton:
% https://github.com/peterhcharlton/pulse-analyse
%
% INPUT: PPG: PPG signal
%       configs: configs struct, see below for details
%       plot_flag: flag whether to plot a summary
% 
% OUTPUT: pts: struct defining the locations, timings and amplitudes of all fiducial points
%         norm_pts: struct defining the locations, timings and amplitudes of all fiducial points using normalised PPG pulses -- only returned if config.do_normalise == 1
%           derivs: struct of PPG derivatives
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
%
% Relevant literature:
% - Charlton, P.H., Celka, P., Farukh, B., Chowienczyk, P. and Alastruey, J., 2018. Assessing mental stress from the photoplethysmogram: a numerical study. Physiological measurement, 39(5), p.054001.
% - Balmer, J., Smith, R., Pretty, C.G., Desaive, T., Shaw, G.M. and Chase, J.G., 2021. Accurate end systole detection in dicrotic notch-less arterial pressure waveforms. Journal of clinical monitoring and computing, 35(1), pp.79-88.
%

narginchk(1, inf);

% What fid points we support
allowed_fid_pt_names = {'a', 'b', 'c', 'd', 'e', 'f', 's', 'dia', 'dic', 'p1pk', ...
    'p2pk', 'p1in', 'p2in', 'W', 'f1', 'f2', 'halfpoint', 'tangent', ...
    'gauss'};
narginchk(1, inf)
if nargin < 2 || isempty(config)
    config = struct();
end
if nargin < 3
    plot_flag  = false;
end
default_config.fid_pt_names = allowed_fid_pt_names;
default_config.gauss_continue_points = 0; % Used in Gaussian fitting
default_config.do_e_for_dic = false; % If to use e as the location of dicrotic notch
default_config.do_filter = 1;
default_config.do_normalise = 1; % Whether to also return the normalised fiducial points
config = func.aux_functions.update_with_default_opts(config, default_config);
fid_pt_names = config.fid_pt_names;
%Do we need to do Gaussian decomposition?
config.do_gauss = any(strcmp(fid_pt_names, 'gauss'));
%% Error checks
if isempty(PPG.ts)
    pts = [];
    derivs = [];
    return
end
%Are there any unsupported fiducial points?
[~,ia] = intersect(fid_pt_names,allowed_fid_pt_names, 'stable');
if length(ia) ~= length(fid_pt_names)
    error('Contains unsupported fiducial points')
end

%sort fid points names so that if tangent is asked it is last (so that f1 is before).
if any(strcmp(fid_pt_names, 'tangent'))
    fid_pt_names(strcmp(fid_pt_names, 'tangent'))= [];
    fid_pt_names{end+1} = 'tangent';
end
%% Filter Signal
if config.do_filter    
    % Butterworth IIR bandpass filter
    [b,a] = butter(8,10/(PPG.fs/2), 'low');
    ts_filt = filtfilt(b,a,PPG.ts);
else
    ts_filt = PPG.ts;
end

%% Calculate Derivatives
derivs = func.pulsew.get_derivs(ts_filt, PPG.fs);
deriv_names = fieldnames(derivs);
%% Identify Fiducial Points, and Calculate Fiducial Point Timings and Amplitudes
num_beats = length(PPG.peaks);

% Set up broadcast storing variables and flags
for fid_pt_no = 1 : length(allowed_fid_pt_names)
    flags.(['do_', allowed_fid_pt_names{fid_pt_no}]) = ismember(allowed_fid_pt_names{fid_pt_no}, fid_pt_names);
    if flags.(['do_', allowed_fid_pt_names{fid_pt_no}])
        store.(allowed_fid_pt_names{fid_pt_no}) = nan(num_beats, 1);
    end
end
starting_index = nan(num_beats,1);
%% Initialise points store
for fid_pt_no = 1 : length(fid_pt_names)
    if strcmp(fid_pt_names{fid_pt_no}, 'gauss')
        continue
    end
    pts.(fid_pt_names{fid_pt_no}).ind = nan(num_beats,1);
    pts.(fid_pt_names{fid_pt_no}).amp = nan(num_beats,1);
    pts.(fid_pt_names{fid_pt_no}).t   = nan(num_beats,1);
    
    if config.do_normalise
        norm_pts.(fid_pt_names{fid_pt_no}).ind = nan(num_beats,1);
        norm_pts.(fid_pt_names{fid_pt_no}).amp = nan(num_beats,1);
        norm_pts.(fid_pt_names{fid_pt_no}).t   = nan(num_beats,1);
    end
end

%% Store half amplitude of PPG pulses
halfpoint_values = PPG.ts(PPG.onsets(1:end-1)) + ...
    0.5*(PPG.ts(PPG.peaks) -...
    PPG.ts(PPG.onsets(1:end-1)));
%% For dicrotic notch detection if using Balmer notch detection
dic_detection_vals.T = PPG.t(PPG.onsets(2:end)) - PPG.t(PPG.onsets(1:end-1));
dic_detection_vals.Tau_func = @(t, t_peak, T)(t - t_peak)/(T - t_peak);
dic_detection_vals.most_recent_t_systole = cell(num_beats, 1);
dic_detection_vals.num_beats_average = 10;
dic_detection_vals.Beta = 5;
%% Run loop
for pulse_no = 1 : num_beats
    % If this pulse has a low SQI then dont detect fiducial points for it
    if PPG.sqi_beat(pulse_no) == 0
        continue
    end
    %% Current PPG pulse that we care about
    % extract data for this pulse wave
    curr_els = PPG.onsets(pulse_no):PPG.onsets(pulse_no+1);
    starting_index(pulse_no) = curr_els(1);
    curr = [];
    curr.t = PPG.t(curr_els);t_start = curr.t(1);
    curr.t = curr.t - t_start;
    curr.ts = ts_filt(curr_els);
    
    % Correct for low frequency baseline drift in a single beat 
    correction_line = linspace(curr.ts(1), curr.ts(end), length(curr.ts));
    curr.ts = curr.ts - correction_line' + curr.ts(1);
    
    curr.ts_orig = curr.ts;    
    
    if config.do_normalise
        curr.ts_norm = curr.ts - min(curr.ts);
        curr.ts_norm = curr.ts_norm/ max(curr.ts_norm);
        curr.fs_norm = length(curr.t);
        curr.t_norm = curr.t / curr.t(end);
        % Get derivatives
        curr.derivs_norm = func.pulsew.get_derivs(curr.ts_norm, curr.fs_norm);
    end
    
    for d_idx = 1:length(deriv_names)
        curr.derivs.(deriv_names{d_idx}) = derivs.(deriv_names{d_idx})(curr_els);
    end
    %% Identify fiducial points    
    %%%%%%%%START%%%%%%%% find f1 and f2
    if flags.do_f1
        store.f1(pulse_no) = 1;
    end
    if flags.do_f2
        store.f2(pulse_no) = length(curr.ts);
    end
    %%%%%%%%END%%%%%%%%
    
    %%%%%%%%START%%%%%%%% find s -- the current peak
    if flags.do_s
        store.s(pulse_no) = PPG.peaks(pulse_no);
        %shift by start of pulse
        store.s(pulse_no) = store.s(pulse_no) - PPG.onsets(pulse_no) +1;
    end
    %%%%%%%%END%%%%%%%%
    
    %%%%%%%%START%%%%%%%% find halfpoint the point that is halfway between
    %%%%%%%%the onset and the peak
    if flags.do_halfpoint
        [~,store.halfpoint(pulse_no)] = min(abs(curr.ts(store.f1(pulse_no) : store.s(pulse_no)) - halfpoint_values(pulse_no)));
    end
    %%%%%%%%END%%%%%%%%
    
    %%%%%%%%START%%%%%%%% find W
    if flags.do_W
        [~, store.W(pulse_no)] = max(curr.derivs.first);
    end
    %%%%%%%%END%%%%%%%%
    
    %%%%%%%%START%%%%%%%% find tangent
    if flags.do_tangent
        %         store.tangent(pulse_no) = ((curr.ts(store.f1(pulse_no))-curr.ts(store.W(pulse_no)))/curr.derivs.first(store.W(pulse_no))) + store.W(pulse_no);
        offset = curr.t(store.W(pulse_no)) - (curr.ts_orig(store.W(pulse_no)) - curr.ts_orig(store.f1(pulse_no)))/curr.derivs.first(store.W(pulse_no));
        store.tangent(pulse_no) = offset*PPG.fs + store.f1(pulse_no);
    end
    %%%%%%%%END%%%%%%%%
    
    %%%%%%%%START%%%%%%%% find a
    if flags.do_a
        pks = func.waveform.find_pks_trs(curr.derivs.second, 'pk'); %local function that points to peaks or troughs of waveform
        rel_pks = pks(find(pks < store.W(pulse_no), 1, 'last')); % a must be before the peak of the VPG
        [~, temp_el] = max(curr.derivs.second(rel_pks));
        temp_a = rel_pks(temp_el);
        
        if isempty(temp_a)
            continue
        end
        store.a(pulse_no) = temp_a;
    end
    %%%%%%%%END%%%%%%%%
    
    
    
    %%%%%%%%START%%%%%%%% find b
    if flags.do_b
        % Find local minima in second derivative
        trs = func.waveform.find_pks_trs(curr.derivs.second, 'tr');
        % define an upper bound as 25% of the duration of the signal
        upper_bound = 0.25*length(curr.ts);
        % find local minima between 'a' and this upper bound
        possible_trs = find(trs > store.a(pulse_no) & curr.derivs.second(trs) < 0 & trs < upper_bound);
        % Identify the lowest of these local minima
        [~, rel_el] = min(curr.derivs.second(trs(possible_trs)));
        possible_trs = possible_trs(rel_el);
        temp_b = trs(possible_trs);
        
        if isempty(temp_b)
            % find local minima after 'a'
            possible_pts = find(trs > store.a(pulse_no) & curr.derivs.second(trs) < 0, 1, 'first');
            temp_b = trs(possible_pts);
            
            if isempty(temp_b)
                continue
            end
        end
        store.b(pulse_no) = temp_b;
    end
       %%%%%%%%END%%%%%%%%
    
    %%%%%%%%START%%%%%%%%
    if flags.do_p1pk || flags.do_p1in
        % find p1 -- early systolic peak
        temp_p1 = func.pulsew.identify_p1(curr, store.b(pulse_no), [0.1, 0.18].* PPG.fs );
    else
        temp_p1 = [];
    end
    %%%%%%%%
    
    %%%%%%%%START%%%%%%%%
    if flags.do_e
        % find e
        temp_e = func.pulsew.identify_e(curr, store.W(pulse_no), store.b(pulse_no));
        %     store.e(pulse_no) = func.pulsew.identify_e(curr, store.W(pulse_no), store.b(pulse_no));
    else
        temp_e = [];
    end
    %%%%%%%%END%%%%%%%%
    
    if isempty(temp_e)
        if flags.do_p2in || flags.do_p2pk
            temp_p2 = deal(nan);
        end
    else
        if flags.do_e
            store.e(pulse_no) = temp_e;
        end
        
        
        %%%%%%%%START%%%%%%%%
        if flags.do_f
            % find f
            lower_bound = store.e(pulse_no);
            trs = func.waveform.find_pks_trs(curr.derivs.second, 'tr');
            possible_els = trs(find(trs >=lower_bound, 1));
            if ~isempty(possible_els)
                store.f(pulse_no) = possible_els(1);
            end
            %         clear possible_els trs lower_bound
        end
        %%%%%%%%END%%%%%%%%
        
        
        %%%%%%%%START%%%%%%%%
        if flags.do_dic 
           if config.do_e_for_dic
                if store.e(pulse_no) < store.s(pulse_no)
                    try
                        if store.dic(pulse_no-1) < store.f2(pulse_no)
                            store.dic(pulse_no) = store.dic(pulse_no-1);
                        else
                            store.dic(pulse_no) = nan;
                        end
                    catch
                        store.dic(pulse_no) = nan;
                    end
                else
                    store.dic(pulse_no) = store.e(pulse_no);
                end
                
            else %Dicrotic notch detection by Balmer et al - used in ICU settings where dicrotic notch is diminished
                t_peak = curr.t(store.s(pulse_no));
                curr.Tau = dic_detection_vals.Tau_func(curr.t, t_peak, dic_detection_vals.T(pulse_no));
                mat_most_recent_t_systole = cell2mat(dic_detection_vals.most_recent_t_systole);
                
                if length(mat_most_recent_t_systole) < dic_detection_vals.num_beats_average
                    t_w_max = 0.45 - 0.1/dic_detection_vals.T(pulse_no);
                else
                    t_w_max = mean(mat_most_recent_t_systole(end-(dic_detection_vals.num_beats_average-1):end));
                end
                Tau_w_max = dic_detection_vals.Tau_func(t_w_max, t_peak, dic_detection_vals.T(pulse_no));
                dic_detection_vals.alpha = (dic_detection_vals.Beta * Tau_w_max  - 2* Tau_w_max +1)/(1 - Tau_w_max );
                dic_detection_vals.alpha = min(max(dic_detection_vals.alpha, 1.5), 4.5);
                curr.w = zeros(size(curr.ts));
                
                loc_start = store.s(pulse_no);
                
                curr.w(loc_start:end) = curr.Tau(loc_start:end).^(dic_detection_vals.alpha - 1)  .* (1 - curr.Tau(loc_start:end)).^(dic_detection_vals.Beta -1);
                
                %0.6T limit -- This is an upper limit of the time to systole
                %as discussed in more depth in the identify_e.m function
                curr.w(round(0.6*(length(curr.ts))):end) = 0;
                
                %             [peak_val,loc_notch]= findpeaks(curr.derivs.second .* curr.w);
                weighted_deriv = curr.derivs.second .* curr.w;
                loc_notch = func.waveform.find_pks_trs(weighted_deriv, 'pks');
                peak_val = weighted_deriv(loc_notch);
                if ~isempty(loc_notch)
                    [~, i_loc_notch] = max(peak_val); 
                    store.dic(pulse_no) = loc_notch(i_loc_notch);
                    
                    %Update most recent t_sys
                    t_w_max = curr.t(store.dic(pulse_no));
                    dic_detection_vals.most_recent_t_systole{pulse_no} = t_w_max;
                end
            end
        end
        %%%%%%%%END%%%%%%%%
        
        
        %%%%%%%%START%%%%%%%%
        if flags.do_dia
            end_buffer = length(curr.ts) - ceil( 0.2 * length(curr.ts));  % in samples
            % find dia
            pks = func.waveform.find_pks_trs(curr.ts, 'pks');
            temp_dia = pks(find(pks > store.dic(pulse_no) & pks < end_buffer, 1));
            %         store.dia(pulse_no) =
            if isempty(temp_dia)
                pks = func.waveform.find_pks_trs(curr.derivs.second, 'pks');
                %                 temp_dia = pks(find(pks > store.e(pulse_no), 1, 'first'));
                temp_dia = pks(find(pks > store.dic(pulse_no), 1, 'first'));
                if ~isempty(temp_dia)
                    store.dia(pulse_no) = temp_dia;
                end
            else
                store.dia(pulse_no) = temp_dia;
            end
            %         clear pks
        end
        %%%%%%%%END%%%%%%%%
        
        
        %%%%%%%%START%%%%%%%%
        if flags.do_c
            % find c
            temp_c = func.pulsew.identify_c(curr, store.b(pulse_no), store.e(pulse_no));
            %         store.c(pulse_no) = func.pulsew.identify_c(curr, store.b(pulse_no), store.e(pulse_no));
        else
            temp_c=[];
        end
        %%%%%%%%END%%%%%%%%
        
        
        if ~isempty(temp_c)
            if flags.do_c
                store.c(pulse_no) = temp_c;
            end
            %%%%%%%%START%%%%%%%%
            if flags.do_d
                % find d - the lowest minimum of the second deriv between "c" and "e"
                trs = func.waveform.find_pks_trs(curr.derivs.second, 'tr');
                possible_trs = find(trs > store.c(pulse_no) & trs < store.e(pulse_no));
                if ~isempty(possible_trs)
                    possible_trs = trs(possible_trs);
                    [~, temp_el] = min(curr.derivs.second(possible_trs));
                    store.d(pulse_no) = possible_trs(temp_el);
                    %                 clear possible_trs
                else
                    % unless there isn't a minimum, in which case it's an inflection, and
                    % "d" is the same as "c"
                    store.d(pulse_no) = store.c(pulse_no);
                end
            end
            %%%%%%%%END%%%%%%%%
            
            
            %%%%%%%%START%%%%%%%%
            if flags.do_p2pk || flags.do_p2in
                % find p2 -- late systolic peak
                temp_p2 = func.pulsew.identify_p2(curr, store.d(pulse_no), temp_p1, store.dic(pulse_no));
            end
            %%%%%%%%END%%%%%%%%
        else
            temp_p2 = nan;
        end
        
    end
    if isempty(temp_p1)
        temp_p1 = nan;
    end
    
    % retain timings of original p1 and p2 estimates
    %Early and late systolic components
    if flags.do_p1in
        store.p1in(pulse_no) = temp_p1;
    end
    if flags.do_p2in
        store.p2in(pulse_no) = temp_p2;
    end
    
    if flags.do_p1pk && flags.do_p2pk
        % make p1 or p2 coincident with the systolic peak
        [~, rel_el] = min(abs(store.s(pulse_no)-[temp_p1,temp_p2]));
        if rel_el == 1
            temp_p1 = store.s(pulse_no);
        else
            temp_p2 = store.s(pulse_no);
        end
        
        if ~isnan(temp_p2) && ~isnan(temp_p1)
            % make sure p2 is at a peak if necessary
            pks = func.waveform.find_pks_trs(curr.ts, 'pk');
            cutoff = mean([temp_p1, temp_p2]);
            possible_pks = find(pks > cutoff & pks < store.e(pulse_no) & curr.ts(pks) > curr.ts(temp_p2));
            if ~isempty(possible_pks)
                [~, temp_el] = max(curr.ts(pks(possible_pks)));
                temp_p2 = pks(possible_pks(temp_el));
            end
            % make sure p1 is at a peak if necessary
            pks = func.waveform.find_pks_trs(curr.ts, 'pk');
            cutoff = mean([temp_p1, temp_p2]);
            possible_pks = find(pks < cutoff & pks > store.W(pulse_no) & curr.ts(pks) > curr.ts(temp_p1));
            if ~isempty(possible_pks)
                [~, temp_el] = max(curr.ts(pks(possible_pks)));
                temp_p1 = pks(possible_pks(temp_el));
            end
        end
        
        % store p1pk and p2pk p1pk == p1 at peak
        store.p1pk(pulse_no) = temp_p1;
        store.p2pk(pulse_no) = temp_p2;
    end
    
    %% Store detected points
    for fid_pt_no = 1 : length(fid_pt_names)
        if strcmp(fid_pt_names{fid_pt_no}, 'gauss')
            continue
        end
        curr_temp_el = store.(fid_pt_names{fid_pt_no})(pulse_no);
        
        %remove nans;
        if isnan(curr_temp_el )
            continue
        end
        
        % - index of fiducial point
        curr_temp_el = curr_temp_el + starting_index(pulse_no)-1;
        pts.(fid_pt_names{fid_pt_no}).ind(pulse_no) = curr_temp_el;
        
        % - amplitude of fiducial point
        if sum(strcmp(fid_pt_names{fid_pt_no}, {'s','dia','dic','p1pk','p1in',...
                'p2pk','p2in','f1','f2', 'halfpoint' }))
            amp = PPG.ts(curr_temp_el);
        elseif strcmp(fid_pt_names{fid_pt_no}, 'tangent')
            %Tangent amp is onset amp
            amp = pts.f1.amp(pulse_no);
        elseif sum(strcmp(fid_pt_names{fid_pt_no}, {'W','ms2'}))
            amp = derivs.first(curr_temp_el);
        elseif   sum(strcmp(fid_pt_names{fid_pt_no}, {'a','b','c','d','e','f'}))
            amp = derivs.second(curr_temp_el);
        else
            error(['Unknown type: ', fid_pt_names{fid_pt_no}])
        end
        pts.(fid_pt_names{fid_pt_no}).amp(pulse_no) = amp;
        
        
        % - timing of fiducial point
        try
            t = PPG.t(curr_temp_el(pulse_no));
        catch
            %This is for tangent that has a fractal index value
            t = interp1(1:length(PPG.t),PPG.t,curr_temp_el);
        end
        pts.(fid_pt_names{fid_pt_no}).t(pulse_no) = t;
        
        
        %%%% NORMALISED pulse
        if config.do_normalise
            curr_temp_el = curr_temp_el  -starting_index(pulse_no) +1;
            norm_pts.(fid_pt_names{fid_pt_no}).ind(pulse_no) = curr_temp_el;
            
            % - amplitude of fiducial point
            if sum(strcmp(fid_pt_names{fid_pt_no}, {'s','dia','dic','p1pk','p1in',...
                    'p2pk','p2in','f1','f2', 'halfpoint' }))
                amp = curr.ts_norm(curr_temp_el);
            elseif strcmp(fid_pt_names{fid_pt_no}, 'tangent')
                %Tangent amp is onset amp
                amp = 0;
            elseif sum(strcmp(fid_pt_names{fid_pt_no}, {'W','ms2'}))
                amp =curr.derivs_norm.first(curr_temp_el);
            elseif   sum(strcmp(fid_pt_names{fid_pt_no}, {'a','b','c','d','e','f'}))
                amp = curr.derivs_norm.second(curr_temp_el);
            else
                error(['Unknown type: ', fid_pt_names{fid_pt_no}])
            end
            norm_pts.(fid_pt_names{fid_pt_no}).amp(pulse_no) = amp;
            
            % - timing of fiducial point
            t = (curr_temp_el - 1)/curr.fs_norm;
            %     t_start = pulse_no -1; % As each beat as been
            norm_pts.(fid_pt_names{fid_pt_no}).t(pulse_no) = t;
        end
        
        
        
    end
end
%% Gaussian profiles
if config.do_gauss
    disp('Now computing Gaussian features -- this may take some time for long recordings')
    gauss_pts = func.pulsew.gaussian_model(ts_filt, PPG.t, PPG.onsets, PPG.sqi_beat, config);
    f = fieldnames(gauss_pts);
    
    for i = 1:length(f)
        if ~config.do_normalise
            pts.(f{i}) = gauss_pts.(f{i});
        else
            norm_pts.(f{i}) = gauss_pts.(f{i});
        end
    end
end
%% Define varargout
varargout{1} = pts;
if config.do_normalise; varargout{end+1} = norm_pts; end
varargout{end+1} = derivs;
%% Now plot
if plot_flag
    % Colours for plotting
    colours = func.aux_functions.define_colours;
    
    markers_store = 'o*^do*^do*^d';
    do_mins = PPG.t(end) > 10*60;
    if do_mins % Plot mins or seconds
        PPG.t =  PPG.t/60;
    end
    
    %% Plot time series with selected fiducial points
    figure('Position', [688   109   794   860])
    
    good_beats = PPG.sqi_beat > 0;
    
    num_rows = 4; num_cols = 1; p_idx = 1;
    
    marker_size = 40;
    face_alpha = 0.3; edge_alpha = 1; edge_color ='k';
    
    
    %PPG plot
    ax(p_idx) = subplot(num_rows, num_cols, p_idx);hold on; p_idx = p_idx+1;
    plot(PPG.t, PPG.ts, 'Color', colours.black, 'LineWidth', 1.2)
    %     end
    m_idx = 1;
    if flags.do_s
        scatter(pts.s.t(good_beats), pts.s.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, ...
            'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color, 'Marker', markers_store(m_idx))
        m_idx = m_idx+1;
    end
    if flags.do_dic
        scatter(pts.dic.t(good_beats), pts.dic.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, ...
            'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color, 'Marker', markers_store(m_idx))
        m_idx = m_idx+1;
    end
    if flags.do_dia
        scatter(pts.dia.t(good_beats), pts.dia.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, ...
            'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color, 'Marker', markers_store(m_idx))
        m_idx = m_idx+1;
    end
    ylabel('PPG')
    %     legend('PPG','Onset', 'Tangent', 'Halfpoint', 'Peak','Dicrotic notch', 'Diastolic peak', 'Orientation', 'horizontal')
    legend('PPG','Peak','Dicrotic notch', 'Diastolic peak', 'Orientation', 'horizontal')
    
    
    
    %VPG plot
    ax(p_idx) = subplot(num_rows, num_cols, p_idx);hold on; p_idx = p_idx+1;
    plot(PPG.t, derivs.first, 'Color', colours.black, 'LineWidth', 1.2)
    scatter(pts.W.t(good_beats), pts.W.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color)
    %    xticks([])
    ylabel('VPG')
    legend('VPG', 'Max slope', 'Orientation', 'horizontal')
    
    %APG plot
    
    ax(p_idx) = subplot(num_rows, num_cols, p_idx);hold on; p_idx = p_idx+1;
    plot(PPG.t, derivs.second, 'Color', colours.black, 'LineWidth', 1.2)
    if flags.do_a; scatter(pts.a.t(good_beats), pts.a.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color); end
    if flags.do_b; scatter(pts.b.t(good_beats), pts.b.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color); end
    if flags.do_c; scatter(pts.c.t(good_beats), pts.c.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color); end
    if flags.do_d; scatter(pts.d.t(good_beats), pts.d.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color); end
    if flags.do_e; scatter(pts.e.t(good_beats), pts.e.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color); end
    if flags.do_f; scatter(pts.f.t(good_beats), pts.f.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color); end
    %    xticks([])
    ylabel('APG')
    legend('APG', 'a', 'b', 'c', 'd', 'e', 'f', 'Orientation', 'horizontal')
    
    %3rd deriv plot
    
    ax(p_idx) = subplot(num_rows, num_cols, p_idx);hold on; p_idx = p_idx+1;
    plot(PPG.t, derivs.third, 'Color', colours.black, 'LineWidth', 1.2)
    if flags.do_p1pk; scatter(pts.p1pk.t(good_beats), pts.p1pk.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color); end
    if flags.do_p2pk; scatter(pts.p2pk.t(good_beats), pts.p2pk.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color); end
    ylabel('JPG')
    legend('JPG', 'p1pk', 'p2pk', 'Orientation', 'horizontal')
    
    if do_mins
        xlabel('Time (mins)');
    else
        xlabel('Time (secs)');
    end
    
    linkaxes(ax, 'x')
    func.plot.tightfig();    
end


end