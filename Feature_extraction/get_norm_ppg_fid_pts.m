% get_ppg_fid_pts - This function locates the fiducial points of the PPG that
% are commonly used in BP estimation.
%Code adapted from P.Charlton:
%https://github.com/peterhcharlton/pulse-analyse
%allowed_fid_pt_names = {'a', 'b', 'c', 'd', 'e', 'f', 's', 'dia', 'dic', 'p1pk', ...
%    'p2pk', 'p1in', 'p2in', 'W', 'f1', 'f2', 'halfpoint', 'tangent' ...
%   'gauss'};
%
%
% LITERATURE:
% "Assessing mental stress from the photoplethysmogram: A numerical study"
% -- P. Charlton


%TO DO improve peak detection by location of fiducial points (for e.g. peak
%cannot come after the dicrotic notch and therefore we must have detected
%the dicrotic peak).
%TO DO -- should we be setting SQI = 0 for beats with missing fiducial
%points?
%TO DO -- there are a few indicies that need to be together
function [pts, norm_ts, derivs ] = get_norm_ppg_fid_pts(PPG,config, plot_flag)
%what fid points do we support?
allowed_fid_pt_names = {'a', 'b', 'c', 'd', 'e', 'f', 's', 'dia', 'dic', 'p1pk', ...
    'p2pk', 'p1in', 'p2in', 'W', 'f1', 'f2', 'halfpoint', 'tangent' ...
    'gauss','skewed_gauss'};
narginchk(1, inf)
if nargin < 2 || isempty(config)
    % setup variables to store fiducial point indices
    fid_pt_names = allowed_fid_pt_names;
else
    if isfield(config, 'fid_pt_names')
        fid_pt_names = config.fid_pt_names;
    else
        fid_pt_names =config;
    end
end
if nargin < 3
    plot_flag  = false;
end

%Are there any unsupported fiducial points?
[~,ia] = intersect(fid_pt_names,allowed_fid_pt_names, 'stable');
if length(ia) ~= length(fid_pt_names)
    error('Contains unsupported fiducial points')
end

%sort fid points names so that if tangent is asked it is last (so that f1
%is before).
if any(strcmp(fid_pt_names, 'tangent'))
    fid_pt_names(strcmp(fid_pt_names, 'tangent'))= [];
    fid_pt_names{end+1} = 'tangent';
end


%Do we need to do Gaussian decomposition?
if any(strcmp(fid_pt_names, 'gauss'))
    do_gauss = 1;
else
    do_gauss = 0;
end
%Do we need to do Skewed Gaussian decomposition?
if any(strcmp(fid_pt_names, 'skewed_gauss'))
    do_skewed_gauss = 1;
else
    do_skewed_gauss = 0;
end
%% Filter Signal
do_filter = 1; % PPG should be filtered from the preprocess_PPG code
if do_filter
    filter.method      = {'IIR' ; 'IIR'};
    filter.type        = {'highpass'; 'lowpass'};
    filter.order       = [8  8];
    filter.fc1         = [0.5, 10]; filter.fc2  = filter.fc1;
    ts_filt =  pe.sigproc.filter.filter(PPG.ts, PPG.fs, filter);
else
    ts_filt = PPG.ts;
end

%% Calculate Derivatives
s_g_filter_len = 9;

%% Identify Fiducial Points, and Calculate Fiducial Point Timings and Amplitudes


%width25 width50 width75

% Set up broadcast variables
flags.do_test = 1;
for fid_pt_no = 1 : length(allowed_fid_pt_names)
    eval(['store_', allowed_fid_pt_names{fid_pt_no} ' = nan(length(PPG.onsets)-1,1);'])
    eval(['flags.do_', allowed_fid_pt_names{fid_pt_no} ' = 0;'])
end
for fid_pt_no = 1 : length(fid_pt_names)
    eval(['flags.do_', fid_pt_names{fid_pt_no} ' = 1;'])
end
starting_indice = nan(length(PPG.onsets)-1,1);
onsets = PPG.onsets;
peaks = PPG.peaks;
sqi_beat = PPG.sqi_beat;

% setup buffer zones
buffer_tol.deb = 0.005; % in secs
buffer_tol.fin = 0.2; % in proportion of beat
buffer_p1 = [0.1, 0.18]; % in secs

halfpoint_values = PPG.ts(PPG.onsets(1:end-1)) + ...
    0.5*(PPG.ts(PPG.peaks) -...
    PPG.ts(PPG.onsets(1:end-1)));


num_beats = length(peaks);


%% initialise points store
for fid_pt_no = 1 : length(fid_pt_names)
    if strcmp(fid_pt_names{fid_pt_no}, 'gauss')
        continue
    end
    eval(['pts.' fid_pt_names{fid_pt_no} '.ind = nan(num_beats,1);'])
    eval(['pts.' fid_pt_names{fid_pt_no} '.amp = nan(num_beats,1);'])
    eval(['pts.' fid_pt_names{fid_pt_no} '.t = nan(num_beats,1);'])
end


norm_ts = cell(num_beats,1);
norm_fs = nan(num_beats,1);
derivs.first = cell(num_beats,1);
derivs.second = cell(num_beats,1);
derivs.third = cell(num_beats,1);
derivs.fourth = cell(num_beats,1);
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
    curr.fs = PPG.fs * curr.t(end);
    curr.t = curr.t / curr.t(end);
    
    dt = 1/curr.fs;
    
    
    curr.ts = ts_filt(curr_els);
    curr.ts = curr.ts - min(curr.ts);
    curr.ts = curr.ts/ max(curr.ts);
    
    
    % Correct for low frequency baseline drift in a single beat
    correction_line = linspace(curr.ts(1), curr.ts(end), length(curr.ts));
    curr.ts = curr.ts - correction_line' + curr.ts(1);
    
    clear correction_line
    
    
    curr.derivs.first = savitzky_golay(curr.ts, 1, s_g_filter_len)./dt;
    curr.derivs.second = savitzky_golay(curr.derivs.first, 1, s_g_filter_len)./dt;
    curr.derivs.third = savitzky_golay(curr.derivs.second, 1, s_g_filter_len)./dt;
    curr.derivs.fourth = savitzky_golay(curr.derivs.third, 1, s_g_filter_len)./dt;
    
    
    norm_ts{pulse_no} = curr.ts;
    derivs.first{pulse_no} = curr.derivs.first;
    derivs.second{pulse_no} = curr.derivs.second;
    derivs.third{pulse_no} = curr.derivs.third;
    derivs.fourth{pulse_no} = curr.derivs.fourth;
    
    
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
    %this is currently different to PC method which looks at the max of the filtered signal -- so will be in a different place to the adusted peak
    if flags.do_s
        store_s(pulse_no) = peaks(pulse_no);
        %shift by start of pulse
        store_s(pulse_no) = store_s(pulse_no) - onsets(pulse_no) +1;
    end
    %%%%%%%%END%%%%%%%%
    
    %%%%%%%%START%%%%%%%% find halfpoint the point that is halfway between
    %%%%%%%%the onset and the peak
    if flags.do_halfpoint
        [~,store_halfpoint(pulse_no)] = min(abs(curr.ts(store_f1(pulse_no) : store_s(pulse_no)) - halfpoint_values(pulse_no)));
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
        p1_buffer = floor(buffer_p1 .* fs);  % in samples
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
            % find dic -- Need to find and implement other algorithms for
            % detecting the dicrotic notch -- This is the algorithm recommended
            % by Elgendi and is used by Mok Ahn et al
            
            %Cant have dicrotic notch before the systlic peak -- need to
            %check though as the systolic peak may be wrong
            if store_e(pulse_no) < store_s(pulse_no)
                try
                    if store_dic(pulse_no-1) < store_f2(pulse_no)
                        store_dic(pulse_no) = store_dic(pulse_no-1);
                    else
                        store_dic(pulse_no) = nan;
                    end
                catch
                    store_dic(pulse_no) = nan;
                end
            else
                store_dic(pulse_no) = store_e(pulse_no);
            end
            
            
        end
        %%%%%%%%
        
        %         % Adjust peak -- the peak cannot come after the dicrotic notch --
        %         % so set to be the highest peak before the dicrotic notch -- this
        %         % will need to be adapted once we discuss where the peak is when 3
        %         % are found.
        %         if flags.do_s
        %             if store_dic(pulse_no) < store_s(pulse_no)
        %                 store_s(pulse_no) =  adjust_peak(curr.ts, store_s(pulse_no),store_dic(pulse_no) );
        %             end
        %         end
        
        %%%%%%%%START%%%%%%%%
        if flags.do_dia
            % find dia
            pks = func.waveform.find_pks_trs(curr.ts, 'pks');
            temp_dia = pks(find(pks > store_dic(pulse_no) & pks < end_buffer, 1));
            %         store_dia(pulse_no) =
            if isempty(temp_dia)
                pks = func.waveform.find_pks_trs(curr.derivs.first, 'pks');
                temp_dia = pks(find(pks > store_e(pulse_no), 1, 'first'));
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
do_normalise = 1;
if do_gauss
    gauss_pts = func.pulsew.gaussian_model(ts_filt, PPG.t, onsets, sqi_beat, do_normalise, do_plot);
    f = fieldnames(gauss_pts);
    for i = 1:length(f)
        pts.(f{i}) = gauss_pts.(f{i});
    end
end

if do_skewed_gauss
    skew_gauss_pts = func.pulsew.skewed_gaussian_model(ts_filt, PPG.t, onsets, sqi_beat, do_normalise, do_plot);
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
% %     include_pulse(~good_loc) = false;
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
    
    
    good_beats = include_pulse;
    
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
    
    bold_figure;
    
    %% Plot 2 - Histogram of all points
    pulse_no = 1;
    func.pulsew.plot_ppg_indices(PPG, pulse_no,   1)
    %% Plot 3
    
    pt_names_plt = fieldnames(pts);
    
    constants_def;
    
    %First work out number of rows and columns
    num_rows = ceil(sqrt(numel(pt_names_plt)));
    num_columns = ceil(sqrt(numel(pt_names_plt)));
    figure
    
    for pt_name_no = 1 : numel(pt_names_plt)
        subplot(num_rows, num_columns, pt_name_no)
        eval(['curr_temp_el = pts.',fid_pt_names{pt_name_no}, '.ind;'])
        %Remove nans
        curr_temp_el(isnan(curr_temp_el)) = [];
        % - amplitude of fiducial point
        if sum(strcmp(fid_pt_names{pt_name_no}, {'s','dia','dic','p1pk','p1in','p2pk','p2in','f1','f2', 'halfpoint', 'tangent'}))
            eval(['amp = pts.',fid_pt_names{pt_name_no}, '.amp;'])
            
        elseif sum(strcmp(fid_pt_names{pt_name_no}, {'a','b','c','d','e','f'}))
            amp = PPG.derivs.second_d(curr_temp_el);
            
        elseif sum(strcmp(fid_pt_names{pt_name_no}, {'W'}))
            amp = PPG.derivs.first_d(curr_temp_el);
        else
            continue
        end
        
        histogram(amp, 'FaceColor', colours.blue)
        title(pt_names_plt{pt_name_no})
    end
    %     if up.options.normalise_pw
    %         suptitle('NOTE: Pulse Wave is normalised')
    %     else
    %         suptitle('Pulse Wave is not normalised')
    %     end
    
    %     func.plot.tightfig();
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 0.6);
    
end


end


function deriv = savitzky_golay(sig, deriv_no, win_size)

%% assign coefficients
% From: https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter#Tables_of_selected_convolution_coefficients
% which are calculated from: A., Gorry (1990). "General least-squares smoothing and differentiation by the convolution (Savitzky?Golay) method". Analytical Chemistry. 62 (6): 570?3. doi:10.1021/ac00205a007.

switch deriv_no
    case 0
        % - smoothing
        switch win_size
            case 5
                coeffs = [-3, 12, 17, 12, -3];
                norm_factor = 35;
            case 7
                coeffs = [-2, 3, 6, 7, 6, 3, -2];
                norm_factor = 21;
            case 9
                coeffs = [-21, 14, 39, 54, 59, 54, 39, 14, -21];
                norm_factor = 231;
            otherwise
                error('Can''t do this window size')
        end
    case 1
        % - first derivative
        switch win_size
            case 5
                coeffs = -2:2;
                norm_factor = 10;
            case 7
                coeffs = -3:3;
                norm_factor = 28;
            case 9
                coeffs = -4:4;
                norm_factor = 60;
            otherwise
                error('Can''t do this window size')
        end
        
    case 2
        % - second derivative
        switch win_size
            case 5
                coeffs = [2,-1,-2,-1,2];
                norm_factor = 7;
            case 7
                coeffs = [5,0,-3,-4,-3,0,5];
                norm_factor = 42;
            case 9
                coeffs = [28,7,-8,-17,-20,-17,-8,7,28];
                norm_factor = 462;
            otherwise
                error('Can''t do this window size')
        end
        
    case 3
        % - third derivative
        switch win_size
            case 5
                coeffs = [-1,2,0,-2,1];
                norm_factor = 2;
            case 7
                coeffs = [-1,1,1,0,-1,-1,1];
                norm_factor = 6;
            case 9
                coeffs = [-14,7,13,9,0,-9,-13,-7,14];
                norm_factor = 198;
            otherwise
                error('Can''t do this window size')
        end
        
    case 4
        % - fourth derivative
        switch win_size
            case 7
                coeffs = [3,-7,1,6,1,-7,3];
                norm_factor = 11;
            case 9
                coeffs = [14,-21,-11,9,18,9,-11,-21,14];
                norm_factor = 143;
            otherwise
                error('Can''t do this window size')
        end
        
    otherwise
        error('Can''t do this order of derivative')
end

if rem(deriv_no, 2) == 1
    coeffs = -1*coeffs;
end

A = [1,0];
filtered_sig = filter(coeffs, A, sig);
s=length(sig);
half_win_size = floor(win_size*0.5);
deriv=[filtered_sig(win_size)*ones(half_win_size,1);filtered_sig(win_size:s);filtered_sig(s)*ones(half_win_size,1)];
deriv = deriv/norm_factor;

end


function peak = adjust_peak(pulse, peak, e)
%find the peaks in pulse
pks = func.waveform.find_pks_trs(pulse, 'pk');

%select the peak that is highest before e
rel_pks  = pks(pks < e);
if isempty(rel_pks)
    return
end

[~, max_el] = max(pulse(rel_pks));
peak = rel_pks(max_el);


end
