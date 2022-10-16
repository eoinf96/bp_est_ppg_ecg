%  Released under the GNU General Public License
%  Copyright (C) 2021  Eoin Finnegan
%  eoin.finnegan@eng.ox.ac.uk
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
function [pts ] = get_norm_ppg_fid_pts(PPG,config, plot_flag)
% This function extracts     
    
    
% get_ppg_fid_pts - This function locates the fiducial points of the PPG that
% are commonly used in BP estimation.
%Code adapted from P.Charlton:
%https://github.com/peterhcharlton/pulse-analyse
%
%what fid points do we support?
allowed_fid_pt_names = {'a', 'b', 'c', 'd', 'e', 'f', 's', 'dia', 'dic', 'p1pk', ...
    'p2pk', 'p1in', 'p2in', 'W', 'f1', 'f2', 'halfpoint', 'tangent', ...
    'gauss','skewed_gauss'};
narginchk(1, inf)
if nargin < 2 || isempty(config)
    % setup variables to store fiducial point indices
    config = struct();
end
if nargin < 3
    plot_flag  = false;
end
default_config.fid_pt_names = allowed_fid_pt_names(1:end-2);
default_config.gauss_continue_points = 0;
default_config.do_e_for_dic = false;
default_config.do_filter = 1;

config = func.aux_functions.update_with_default_opts(config, default_config);
fid_pt_names = config.fid_pt_names;
%Do we need to do Gaussian decomposition?
config.do_gauss = any(strcmp(fid_pt_names, 'gauss'));
%Do we need to do Skewed Gaussian decomposition?
config.do_skewed_gauss = any(strcmp(fid_pt_names, 'skewed_gauss'));
%% Error check
if isempty(PPG.ts)
    pts = [];
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
%% Identify Fiducial Points, and Calculate Fiducial Point Timings and Amplitudes
num_beats = length(PPG.peaks);

% Set up broadcast storing variables and flags
for fid_pt_no = 1 : length(allowed_fid_pt_names)
    flags.(['do_', allowed_fid_pt_names{fid_pt_no}]) = ismember(allowed_fid_pt_names{fid_pt_no}, fid_pt_names);
    if flags.(['do_', allowed_fid_pt_names{fid_pt_no}])
        store.(allowed_fid_pt_names{fid_pt_no}) = nan(num_beats, 1);
    end
end
starting_indice = nan(length(PPG.onsets)-1,1);

% setup buffer zones
buffer_tol.deb = 0.005; % in secs
buffer_tol.fin = 0.2; % in proportion of beat % Used for what end proportion of beat wont contain diastolic peak
buffer_p1 = [0.1, 0.18]; % in secs
%% For dicrotic notch detection if using Balmer notch detection
dic_detection_vals.T = PPG.t(PPG.onsets(2:end)) - PPG.t(PPG.onsets(1:end-1));
dic_detection_vals.Tau_func = @(t, t_peak, T)(t - t_peak)/(T - t_peak);
dic_detection_vals.most_recent_t_systole = cell(num_beats, 1);
dic_detection_vals.num_beats_average = 10;
dic_detection_vals.Beta = 5;
%% cycle through each pulse wave
% parfor pulse_no = 1 : length(PPG.onsets)-1
for pulse_no = 1 : num_beats
    
    
    % If this pulse has a low SQI then dont detect fiducial points for it
    if sqi_beat(pulse_no) == 0
        continue
    end
    
    %%
    % extract data for this pulse wave
    curr_els = onsets(pulse_no):onsets(pulse_no+1);
    starting_indice(pulse_no) = curr_els(1);
    
    curr = [];
    curr.t = PPG.t(curr_els);
    curr.t  = curr.t - curr.t(1);
    curr.fs = length(curr.t);
    curr.t = curr.t / curr.t(end);
    
    curr.ts = ts_filt(curr_els);
    curr.ts = curr.ts - min(curr.ts);
    curr.ts = curr.ts/ max(curr.ts);
    
    
    % Correct for low frequency baseline drift in a single beat
    correction_line = linspace(curr.ts(1), curr.ts(end), length(curr.ts));
    curr.ts = curr.ts - correction_line' + curr.ts(1);    
    % Get derivatives
    curr.derivs = func.pulsew.get_derivs(curr.ts, curr.fs);

    
    %% Identify fiducial points
    
    %%%%%%%%START%%%%%%%% find buffers
    end_buffer = length(curr.ts) - ceil(buffer_tol.fin * length(curr.ts));  % in samples
    %%%%%%%%END%%%%%%%%
    
    %%%%%%%%START%%%%%%%% find f1 and f2
    if flags.do_f1 ==1
        store_f1(pulse_no) = 1;
    end
    if flags.do_f2
        store_f2(pulse_no) = length(curr.ts);
    end
    %%%%%%%%END%%%%%%%%
    
    %%%%%%%%START%%%%%%%% find s -- the current peak
    if flags.do_s
        store_s(pulse_no) = peaks(pulse_no);
        %shift by start of pulse
        store_s(pulse_no) = store_s(pulse_no) - onsets(pulse_no) +1;
    end
    %%%%%%%%END%%%%%%%%
    
    %%%%%%%%START%%%%%%%% find halfpoint the point that is halfway between
    %%%%%%%%the onset and the peak
    if flags.do_halfpoint
        [~,store_halfpoint(pulse_no)] = min(abs(curr.ts(store_f1(pulse_no) : store_s(pulse_no)) - 0.5));
        
    end
    %%%%%%%%END%%%%%%%%
    
    %%%%%%%%START%%%%%%%% find W
    if flags.do_W
        [~, store_W(pulse_no)] = max(curr.derivs.first);
    end
    %%%%%%%%END%%%%%%%%
    
    %%%%%%%%START%%%%%%%% find tangent
    if flags.do_tangent
        %         store_tangent(pulse_no) = ((curr.ts(store_f1(pulse_no))-curr.ts(store_W(pulse_no)))/curr.derivs.first(store_W(pulse_no))) + store_W(pulse_no);
        offset = ( curr.ts(store_f1(pulse_no)) +  curr.derivs.first(store_W(pulse_no)) * curr.t(store_W(pulse_no)) - curr.ts(store_W(pulse_no)))/curr.derivs.first(store_W(pulse_no));
        store_tangent(pulse_no) = offset*curr.fs + store_f1(pulse_no);
    end
    %%%%%%%%END%%%%%%%%
    
    %%%%%%%%START%%%%%%%% find a
    if flags.do_a
        pks = func.waveform.find_pks_trs(curr.derivs.second, 'pk'); %local function that points to peaks or troughs of waveform
        rel_pks = pks(find(pks < store_W(pulse_no), 1, 'last')); % a must be before the peak of the VPG
        [~, temp_el] = max(curr.derivs.second(rel_pks));
        temp_a = rel_pks(temp_el);
        %     clear temp_el pks rel_pks
        
        if isempty(temp_a)
            continue
        end
        store_a(pulse_no) = temp_a;
    end
    %%%%%%%%END%%%%%%%%
    
    
    
    %%%%%%%%START%%%%%%%% find b
    if flags.do_b
        % Find local minima in second derivative
        trs = func.waveform.find_pks_trs(curr.derivs.second, 'tr');
        % define an upper bound as 25% of the duration of the signal
        upper_bound = 0.25*length(curr.ts);
        % find local minima between 'a' and this upper bound
        possible_trs = find(trs > store_a(pulse_no) & curr.derivs.second(trs) < 0 & trs < upper_bound);
        % Identify the lowest of these local minima
        [~, rel_el] = min(curr.derivs.second(trs(possible_trs)));
        possible_trs = possible_trs(rel_el);
        temp_b = trs(possible_trs);
        
        if isempty(temp_b)
            % find local minima after 'a'
            possible_pts = find(trs > store_a(pulse_no) & curr.derivs.second(trs) < 0, 1, 'first');
            temp_b = trs(possible_pts);
            
            if isempty(temp_b)
                continue
            end
        end
        store_b(pulse_no) = temp_b;
    end
    
    %  clear upper_bound  rel_el trs temp
    %%%%%%%%END%%%%%%%%
    
    %%%%%%%%START%%%%%%%%
    if flags.do_p1pk || flags.do_p1in
        % find p1 -- early systolic peak
        p1_buffer = floor(buffer_p1 .* PPG.fs);  % in samples
        temp_p1 = func.pulsew.identify_p1(curr, store_b(pulse_no), p1_buffer);
    else
        temp_p1 = [];
    end
    %%%%%%%%
    
    %%%%%%%%START%%%%%%%%
    if flags.do_e
        % find e
        temp_e = func.pulsew.identify_e(curr, store_W(pulse_no), store_b(pulse_no));
        %     store_e(pulse_no) = func.pulsew.identify_e(curr, store_W(pulse_no), store_b(pulse_no));
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
            store_e(pulse_no) = temp_e;
        end
        
        
        %%%%%%%%START%%%%%%%%
        if flags.do_f
            % find f
            lower_bound = store_e(pulse_no);
            trs = func.waveform.find_pks_trs(curr.derivs.second, 'tr');
            possible_els = trs(find(trs >=lower_bound, 1));
            if ~isempty(possible_els)
                store_f(pulse_no) = possible_els(1);
            end
            %         clear possible_els trs lower_bound
        end
        %%%%%%%%END%%%%%%%%
        
        
        %%%%%%%%START%%%%%%%%
        if flags.do_dic
            % Dicrotic notch detection by Balmer et al
            %Can maybe toggle based on if MIMIC or not
            
            % find dic -- Need to find and implement other algorithms for
            % detecting the dicrotic notch -- This is the algorithm recommended
            % by Elgendi and is used by Mok Ahn et al
            
            %Cant have dicrotic notch before the sysotlic peak -- need to
            %check though as the systolic peak may be wrong
%             if store_e(pulse_no) < store_s(pulse_no)
%                 try
%                     if store_dic(pulse_no-1) < store_f2(pulse_no)
%                         store_dic(pulse_no) = store_dic(pulse_no-1);
%                     else
%                         store_dic(pulse_no) = nan;
%                     end
%                 catch
%                     store_dic(pulse_no) = nan;
%                 end
%             else
%                 store_dic(pulse_no) = store_e(pulse_no);
%             end
           
            t_peak = curr.t(store_s(pulse_no));
            try
            curr.Tau = dic_detection_vals.Tau_func(curr.t, t_peak, dic_detection_vals.T(pulse_no));
            catch
                a=1;
            end
            mat_most_recent_t_systole = cell2mat(dic_detection_vals.most_recent_t_systole);
            if length(mat_most_recent_t_systole) < dic_detection_vals.num_beats_average
                t_w_max = 0.45 - 0.1/dic_detection_vals.T(pulse_no);
            else
                % Need to sort this out -- what if t_sys is wrong?
                t_w_max = mean(mat_most_recent_t_systole(end-(dic_detection_vals.num_beats_average-1):end));
            end
            Tau_w_max = dic_detection_vals.Tau_func(t_w_max, t_peak, dic_detection_vals.T(pulse_no));
            dic_detection_vals.alpha = (dic_detection_vals.Beta * Tau_w_max  - 2* Tau_w_max +1)/(1 - Tau_w_max );
            dic_detection_vals.alpha = min(max(dic_detection_vals.alpha, 1.5), 4.5);
            curr.w = zeros(size(curr.ts));
        
            loc_start = store_s(pulse_no);

            curr.w(loc_start:end) = curr.Tau(loc_start:end).^(dic_detection_vals.alpha - 1)  .* (1 - curr.Tau(loc_start:end)).^(dic_detection_vals.Beta -1);

            %0.6T limit -- This is an upper limit of the time to systole
            %as discussed in more depth in the identify_e.m function
            curr.w(round(0.6*(length(curr.ts))):end) = 0;

%             [peak_val,loc_notch]= findpeaks(curr.derivs.second .* curr.w);
            weighted_deriv = curr.derivs.second .* curr.w;
            loc_notch = func.waveform.find_pks_trs(weighted_deriv, 'pks');
            peak_val = weighted_deriv(loc_notch);
            if ~isempty(loc_notch)
                [~, i_loc_notch] = max(peak_val); %%% he maximum peak??
                store_dic(pulse_no) = loc_notch(i_loc_notch);
                
                %Update most recent t_sys
                t_w_max = curr.t(store_dic(pulse_no));
                dic_detection_vals.most_recent_t_systole{pulse_no} = t_w_max;
            end
            %Can also move the notch to the point of VPG closest to 0. Not
            %done here but check Balmer_notch_detection for implementation.
        end
        
        
        %%%%%%%%

        
        %%%%%%%%START%%%%%%%%
        if flags.do_dia
            % find dia
            pks = func.waveform.find_pks_trs(curr.ts, 'pks');
            temp_dia = pks(find(pks > store_dic(pulse_no) & pks < end_buffer, 1));
            %         store_dia(pulse_no) =
            if isempty(temp_dia)
                pks = func.waveform.find_pks_trs(curr.derivs.first, 'pks');
                temp_dia = pks(find(pks > store_dic(pulse_no), 1, 'first'));
                if ~isempty(temp_dia)
                    store_dia(pulse_no) = temp_dia;
                end
            else
                store_dia(pulse_no) = temp_dia;
            end
            %         clear pks
        end
        %%%%%%%%END%%%%%%%%
        
        
        %%%%%%%%START%%%%%%%%
        if flags.do_c
            % find c
            temp_c = func.pulsew.identify_c(curr, store_b(pulse_no), store_e(pulse_no));
            %         store_c(pulse_no) = func.pulsew.identify_c(curr, store_b(pulse_no), store_e(pulse_no));
        else
            temp_c=[];
        end
        %%%%%%%%END%%%%%%%%
        
        
        if ~isempty(temp_c)
            if flags.do_c
                store_c(pulse_no) = temp_c;
            end
            %%%%%%%%START%%%%%%%%
            if flags.do_d
                % find d - the lowest minimum of the second deriv between "c" and "e"
                trs = func.waveform.find_pks_trs(curr.derivs.second, 'tr');
                possible_trs = find(trs > store_c(pulse_no) & trs < store_e(pulse_no));
                if ~isempty(possible_trs)
                    possible_trs = trs(possible_trs);
                    [~, temp_el] = min(curr.derivs.second(possible_trs));
                    store_d(pulse_no) = possible_trs(temp_el);
                    %                 clear possible_trs
                else
                    % unless there isn't a minimum, in which case it's an inflection, and
                    % "d" is the same as "c"
                    store_d(pulse_no) = store_c(pulse_no);
                end
            end
            %%%%%%%%END%%%%%%%%
            
            
            %%%%%%%%START%%%%%%%%
            if flags.do_p2pk || flags.do_p2in
                % find p2 -- late systolic peak
                temp_p2 = func.pulsew.identify_p2(curr, store_d(pulse_no), temp_p1, store_dic(pulse_no));
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
        store_p1in(pulse_no) = temp_p1;
    end
    if flags.do_p2in
        store_p2in(pulse_no) = temp_p2;
    end
    
    if flags.do_p1pk && flags.do_p2pk
        % make p1 or p2 coincident with the systolic peak
        [~, rel_el] = min(abs(store_s(pulse_no)-[temp_p1,temp_p2]));
        if rel_el == 1
            temp_p1 = store_s(pulse_no);
        else
            temp_p2 = store_s(pulse_no);
        end
        
        if ~isnan(temp_p2) && ~isnan(temp_p1)
            % make sure p2 is at a peak if necessary
            pks = func.waveform.find_pks_trs(curr.ts, 'pk');
            cutoff = mean([temp_p1, temp_p2]);
            possible_pks = find(pks > cutoff & pks < store_e(pulse_no) & curr.ts(pks) > curr.ts(temp_p2));
            if ~isempty(possible_pks)
                [~, temp_el] = max(curr.ts(pks(possible_pks)));
                temp_p2 = pks(possible_pks(temp_el));
            end
            % make sure p1 is at a peak if necessary
            pks = func.waveform.find_pks_trs(curr.ts, 'pk');
            cutoff = mean([temp_p1, temp_p2]);
            possible_pks = find(pks < cutoff & pks > store_W(pulse_no) & curr.ts(pks) > curr.ts(temp_p1));
            if ~isempty(possible_pks)
                [~, temp_el] = max(curr.ts(pks(possible_pks)));
                temp_p1 = pks(possible_pks(temp_el));
            end
        end
        
        % store p1pk and p2pk p1pk == p1 at peak
        store_p1pk(pulse_no) = temp_p1;
        store_p2pk(pulse_no) = temp_p2;
        
    end
    
    %% Save each fid point on the normalised beat
    for fid_pt_no = 1 : length(fid_pt_names)
        
        
        eval(['curr_temp_el = store_' fid_pt_names{fid_pt_no} '(pulse_no);']);
        
        if isnan(curr_temp_el )
            continue
        end
        % - index of fiducial point
        %     curr_temp_el = curr_temp_el + starting_indice(pulse_no)-1;
        %     eval(['pts.' fid_pt_names{fid_pt_no} '.ind(pulse_no) = curr_temp_el + starting_indice(pulse_no)-1;'])
        eval(['pts.' fid_pt_names{fid_pt_no} '.ind(pulse_no) = curr_temp_el;'])
        
        
        % - amplitude of fiducial point
        if sum(strcmp(fid_pt_names{fid_pt_no}, {'s','dia','dic','p1pk','p1in',...
                'p2pk','p2in','f1','f2', 'halfpoint' }))
            amp = curr.ts(curr_temp_el);
        elseif strcmp(fid_pt_names{fid_pt_no}, 'tangent')
            %Tangent amp is onset amp
            amp = 0;
        elseif sum(strcmp(fid_pt_names{fid_pt_no}, {'W','ms2'}))
            amp =curr.derivs.first(curr_temp_el);
        elseif   sum(strcmp(fid_pt_names{fid_pt_no}, {'a','b','c','d','e','f'}))
            amp = curr.derivs.second(curr_temp_el);
        else
            error(['Unknown type: ', fid_pt_names{fid_pt_no}])
        end
        eval(['pts.' fid_pt_names{fid_pt_no} '.amp(pulse_no) = amp;'])
        
        
        
        % - timing of fiducial point
        t = (curr_temp_el - 1)/curr.fs;
        %     t_start = pulse_no -1; % As each beat as been
        eval(['pts.' fid_pt_names{fid_pt_no} '.t(pulse_no) = t;'])
        
        %
        
    end
    
end



%% Add Gaussian profiles

do_plot = 0;
gaussparams.do_normalise = 1;
if ~isfield(config, 'gauss_continue_points')
    gaussparams.continue_points = 0;
else
    gaussparams.continue_points = config.gauss_continue_points;
end
if do_gauss
    disp('Now computing Gaussian features -- this may take some time for long recordings')
    gauss_pts = func.pulsew.gaussian_model(ts_filt, PPG.t, onsets, sqi_beat, gaussparams, do_plot);
    f = fieldnames(gauss_pts);
    for i = 1:length(f)
        pts.(f{i}) = gauss_pts.(f{i});
    end
end

if do_skewed_gauss
    skew_gauss_pts = func.pulsew.skewed_gaussian_model(ts_filt, PPG.t, onsets, sqi_beat, gaussparams, do_plot);
    f = fieldnames(skew_gauss_pts);
    for i = 1:length(f)
        pts.(f{i}) = skew_gauss_pts.(f{i});
    end
end
%%
% %we only want to keep pulses that never have a nan for an index
% include_pulse = true(length(onsets)-1,1);
%
% % store points
% for fid_pt_no = 1 : length(fid_pt_names)
%
%     if strcmp(fid_pt_names{fid_pt_no}, 'gauss')
%         continue
%     end
%
%
%     eval(['curr_temp_el = store_' fid_pt_names{fid_pt_no} ';']);
%
%
%     %remove nans;
%     good_loc = ~isnan(curr_temp_el);
%     if any(~good_loc)
%        a=1;
%     end
%     include_pulse(~good_loc) = false;
%     curr_temp_el=curr_temp_el(good_loc);
%
%
%     % - index of fiducial point
%     curr_temp_el = curr_temp_el + starting_indice(good_loc)-1;
%     eval(['pts.' fid_pt_names{fid_pt_no} '.ind(good_loc) = curr_temp_el;'])
%     eval(['clear store_' fid_pt_names{fid_pt_no}])
%
%
%     % - amplitude of fiducial point
%     if sum(strcmp(fid_pt_names{fid_pt_no}, {'s','dia','dic','p1pk','p1in',...
%             'p2pk','p2in','f1','f2', 'halfpoint' }))
%         amp = PPG.ts(curr_temp_el);
%     elseif strcmp(fid_pt_names{fid_pt_no}, 'tangent')
%         %Tangent amp is onset amp
%         amp = pts.f1.amp(good_loc);
%     elseif sum(strcmp(fid_pt_names{fid_pt_no}, {'W','ms2'}))
%         amp = PPG.derivs.first_d(curr_temp_el);
%     elseif   sum(strcmp(fid_pt_names{fid_pt_no}, {'a','b','c','d','e','f'}))
%         amp = PPG.derivs.second_d(curr_temp_el);
%     else
%         error(['Unknown type: ', fid_pt_names{fid_pt_no}])
%     end
%     eval(['pts.' fid_pt_names{fid_pt_no} '.amp(good_loc) = amp;'])
%
%
%
%     % - timing of fiducial point
%     t = (curr_temp_el-1)/PPG.fs;
%     eval(['pts.' fid_pt_names{fid_pt_no} '.t(good_loc) = t;'])
%
%     %
%     clear amp_norm sig_ind amp t
% end
%
%
% clear starting_indice curr_temp_el
%
% %% Ensure only pts are provided for those pulses with all pts available
% for pt_name_no = 1 : length(fid_pt_names)
%     if strcmp(fid_pt_names{pt_name_no}, 'gauss') || strcmp(fid_pt_names{pt_name_no}, 'f1') || strcmp(fid_pt_names{pt_name_no}, 'f2')
%         continue
%     end
%     eval(['pts.' fid_pt_names{pt_name_no} '.ind(~include_pulse) = nan;']);
%     eval(['pts.' fid_pt_names{pt_name_no} '.amp(~include_pulse) = nan;']);
%     eval(['pts.' fid_pt_names{pt_name_no} '.t(~include_pulse) = nan;']);
% end

%Set SQI of PPG to 0 where all fiducial points have not been able to be
%located
% PPG.sqi_beat(~include_pulse) = 0;

PPG.fid_pts = pts;
%% Now plot


if plot_flag
    %% Plot 1 - time series with selected fiducial points
    figure
    % Colours for plotting
    colours = constants_def('Colors');
    
    
    good_beats = PPG.sqi_beat > 0.8;
    time = PPG.t;
    num_subplots = 4;
    
    marker_size = 40;
    face_alpha = 0.3; edge_alpha = 1; edge_color ='k';
    
    
    %PPG plot
    
    ax(1) = subplot(num_subplots,1,1);
    plot(PPG.t, PPG.ts, 'Color', colours.blue, 'LineWidth', 1.2)
    hold on
    if flags.do_f1
        scatter(pts.f1.t(good_beats), pts.f1.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color)
    end
    if flags.do_tangent
        scatter(pts.tangent.t(good_beats), pts.tangent.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color)
    end
    if flags.do_halfpoint
        scatter(pts.halfpoint.t(good_beats), pts.halfpoint.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color)
    end
    if flags.do_s
        scatter(pts.s.t(good_beats), pts.s.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color)
    end
    if flags.do_dic
        scatter(pts.dic.t(good_beats), pts.dic.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color)
    end
    if flags.do_dia
        scatter(pts.dia.t(good_beats), pts.dia.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color)
    end
    %    xticks([])
    ylabel('PPG')
    legend('PPG','Onset', 'Tangent', 'Halfpoint', 'Peak','Dicrotic notch', 'Diastolic peak', 'Orientation', 'horizontal')
    grid on
    
    
    %VPG plot
    ax(2) = subplot(num_subplots,1,2);
    plot(time, PPG.derivs.first_d, 'Color', colours.blue, 'LineWidth', 1.2)
    hold on
    scatter(pts.W.t(good_beats), pts.W.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color)
    %    xticks([])
    ylabel('VPG')
    legend('VPG', 'Max slope', 'Orientation', 'horizontal')
    grid on
    
    %APG plot
    
    ax(3) = subplot(num_subplots,1,3);
    plot(time, PPG.derivs.second_d, 'Color', colours.blue, 'LineWidth', 1.2)
    hold on
    if flags.do_a; scatter(pts.a.t(good_beats), pts.a.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color); end
    if flags.do_b; scatter(pts.b.t(good_beats), pts.b.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color); end
    if flags.do_c; scatter(pts.c.t(good_beats), pts.c.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color); end
    if flags.do_d; scatter(pts.d.t(good_beats), pts.d.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color); end
    if flags.do_e; scatter(pts.e.t(good_beats), pts.e.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color); end
    if flags.do_f; scatter(pts.f.t(good_beats), pts.f.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color); end
    %    xticks([])
    ylabel('APG')
    legend('APG', 'a', 'b', 'c', 'd', 'e', 'f', 'Orientation', 'horizontal')
    grid on
    
    %3rd deriv plot
    
    ax(4) = subplot(num_subplots,1,4);
    plot(time, PPG.derivs.third_d, 'Color', colours.blue, 'LineWidth', 1.2)
    hold on
    if flags.do_p1pk; scatter(pts.p1pk.t(good_beats), pts.p1pk.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color); end
    if flags.do_p2pk; scatter(pts.p2pk.t(good_beats), pts.p2pk.amp(good_beats), marker_size, 'filled', 'MarkerFaceAlpha', face_alpha, 'MarkerEdgeAlpha', edge_alpha, 'MarkerEdgeColor', edge_color); end
    ylabel('JPG')
    legend('JPG', 'p1pk', 'p2pk', 'Orientation', 'horizontal')
    grid on
    
    xlabel('Time (s)')
    
    linkaxes(ax, 'x')
    
    
%     %% Plot 2 - Histogram of all points
%     pulse_no = 1;
%     func.pulsew.plot_ppg_indices(PPG, pulse_no,   1)
%     %% Plot 3
%     
%     pt_names_plt = fieldnames(pts);
%     
%     constants_def;
%     
%     %First work out number of rows and columns
%     num_rows = ceil(sqrt(numel(pt_names_plt)));
%     num_columns = ceil(sqrt(numel(pt_names_plt)));
%     figure
%     
%     for pt_name_no = 1 : numel(pt_names_plt)
%         subplot(num_rows, num_columns, pt_name_no)
%         eval(['curr_temp_el = pts.',fid_pt_names{pt_name_no}, '.ind;'])
%         %Remove nans
%         curr_temp_el(isnan(curr_temp_el)) = [];
%         % - amplitude of fiducial point
%         if sum(strcmp(fid_pt_names{pt_name_no}, {'s','dia','dic','p1pk','p1in','p2pk','p2in','f1','f2', 'halfpoint', 'tangent'}))
%             eval(['amp = pts.',fid_pt_names{pt_name_no}, '.amp;'])
%             
%         elseif sum(strcmp(fid_pt_names{pt_name_no}, {'a','b','c','d','e','f'}))
%             amp = PPG.derivs.second_d(curr_temp_el);
%             
%         elseif sum(strcmp(fid_pt_names{pt_name_no}, {'W'}))
%             amp = PPG.derivs.first_d(curr_temp_el);
%         else
%             continue
%         end
%         
%         histogram(amp, 'FaceColor', colours.blue)
%         title(pt_names_plt{pt_name_no})
%     end
%     %     if up.options.normalise_pw
%     %         suptitle('NOTE: Pulse Wave is normalised')
%     %     else
%     %         suptitle('Pulse Wave is not normalised')
%     %     end
%     
%     %     func.plot.tightfig();
%     set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 0.6);
%     
end


end
