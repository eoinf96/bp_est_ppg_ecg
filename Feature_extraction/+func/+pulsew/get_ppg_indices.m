% get_ppg_indices - This function computes the inicies from the PPG that
% are commonly used in BP estimation.
%Code adapted from P.Charlton: https://github.com/peterhcharlton/RRest

% Inputs : PPG -- This struct must contain all relevent fiducial points
% Outputs : pw_inds -- a struct of pulse wave indices computed for beat located
% in the PPG including the median value computed from beats of high signal
% quality
function [pw_inds] = get_ppg_indices(PPG,  plot_flag,database, do_filter, sqi_threshold)
narginchk(1, inf);

if nargin < 2
    plot_flag = false;
end
if nargin < 3 || isempty(database)
    ht = nan;
else
    configs = constants_def(database);
    if strcmp(database, 'MOLLIE')
        volunteer_idx = PPG.record_name(2:3);
        volunteer_idx = str2double(volunteer_idx);
        ht = configs.demographics.height(volunteer_idx);
    elseif strcmpi(database, 'Ambulatory')
        volunteer_id = PPG.session_name;
        ht = configs.demographics.height(strcmp(configs.root_session_names, volunteer_id(1:end-1)));
    end
end
if nargin < 4
    do_filter = 1;
end
if nargin < 5
    sqi_threshold = 0.8;
end

if ~isfield(PPG, 'fid_pts')
    PPG.fid_pts= func.pulsew.get_ppg_fid_pts(PPG);
end
fid_pts = PPG.fid_pts;

if isempty(fid_pts)
    pw_inds = [];
    return
end

do_gauss = 0;
if isfield(PPG.fid_pts, 'g1')
    do_gauss = 1;
end
verbose_flag = 0;

num_beats = length(PPG.peaks);


%% Filter Signal
if do_filter
    % Butterworth IIR bandpass filter
    [b,a] = butter(8,10/(PPG.fs/2), 'low');
    ts_filt = filtfilt(b,a,PPG.ts);
else
    ts_filt = PPG.ts;
end

%% Check if derivatives exist and if not then create
if ~isfield(PPG, 'derivs')
s_g_filter_len = 7;
dt = 1/PPG.fs;
PPG.derivs.first_d  = func.waveform.savitzky_golay_deriv(ts_filt, 1, s_g_filter_len)./dt;
PPG.derivs.second_d = func.waveform.savitzky_golay_deriv(ts_filt, 2, s_g_filter_len)./dt./dt;
PPG.derivs.third_d  = func.waveform.savitzky_golay_deriv(ts_filt, 3, s_g_filter_len)./dt./dt./dt;
PPG.derivs.fourth_d = func.waveform.savitzky_golay_deriv(ts_filt, 4, s_g_filter_len)./dt./dt./dt./dt;
end
%% Timings
if verbose_flag
    fprintf('Cannot get ...\n')
end

% - Time between systolic and diastolic peaks (secs)
try
    pw_inds.delta_t = fid_pts.dia.t-fid_pts.s.t;
catch
    if verbose_flag
        fprintf('               - Delta T \n')
    end
end

% - Crest time (secs)
try
    pw_inds.CT = fid_pts.s.t-fid_pts.f1.t;
catch
    if verbose_flag
        fprintf('               - Crest Time \n')
    end
end
% - Crest time divided by height (s/m)
% pw_inds.CT_div_ht = pw_inds.CT./ht;

% - Stiffness index (m/s)
% pw_inds.SI = ht./pw_inds.delta_t;

% - Duration of systole
try
    pw_inds.t_sys = fid_pts.dic.t-fid_pts.f1.t;
catch
    if verbose_flag
        fprintf('               - Duration of systole \n')
    end
end

% - Duration of diastole
try
    pw_inds.t_dia = fid_pts.f2.t-fid_pts.dic.t;
catch
    if verbose_flag
        fprintf('               - Duration of diastole \n')
    end
end

% - Ratio of systolic to diastolic durations
try
    pw_inds.t_ratio = pw_inds.t_sys./pw_inds.t_dia;
catch
    if verbose_flag
        fprintf('               - Ratio of systolic to diastolic durations \n')
    end
end


% - Time from p1 to diastolic peak
% pw_inds.t_p1in_dia = fid_pts.t.dia-fid_pts.t.p1in;
% pw_inds.t_p1pk_dia = fid_pts.t.dia-fid_pts.t.p1pk;

% - Time from p2 to diastolic peak
% pw_inds.t_p2in_dia = fid_pts.t.dia-fid_pts.t.p2in;
% pw_inds.t_p2pk_dia = fid_pts.t.dia-fid_pts.t.p2pk;

% - Time from 'b' to 'c'
% pw_inds.t_b_c = fid_pts.t.c-fid_pts.t.b;

% - Time from 'b' to 'd'
% pw_inds.t_b_d = fid_pts.t.d-fid_pts.t.b;

%% Amplitudes

% - Pulse amplitude
% pw_inds.pulse_amp = fid_pts.amp.s-fid_pts.amp.f1;
% pw_inds.pulse_amp_p1 = fid_pts.amp.p1in-fid_pts.amp.f1;
% pw_inds.pulse_amp_p2 = fid_pts.amp.p2in-fid_pts.amp.f1;
try
    pulse_amp = fid_pts.s.amp-fid_pts.f1.amp;
catch
    if verbose_flag
        fprintf('               - Pulse Amplitude \n')
    end
end

% - Augmentation Pressure
% try
%     pw_inds.AP = fid_pts.p2pk.amp-fid_pts.p1in.amp;
% catch
%     if verbose_flag
%         fprintf('               - Augmentation Pressure \n')
%     end
% end


% % - Agumentation Index Takazawa
% try
%     pw_inds.AI = PPG.fid_pts.dia.amp./PPG.fid_pts.s.amp;
% catch
%     if verbose_flag
%         fprintf('               - Augmentation Index \n')
%     end
% end

% - Dicrotic notch amplitude
try
    pw_inds.dic_amp = fid_pts.dic.amp;
catch
    if verbose_flag
        fprintf('               - Dicrotic notch amplitude \n')
    end
end

% - Reflection Index (calculated using systolic peak)
try
    pw_inds.RI = (fid_pts.dia.amp - fid_pts.f1.amp);
catch
    if verbose_flag
        fprintf('               - Reflection index \n')
    end
end

% - Slope Transit Time (from Addison et al)
try
    pw_inds.STT = (fid_pts.s.t - fid_pts.f1.t)./(fid_pts.s.amp - fid_pts.f1.amp);
catch
    if verbose_flag
        fprintf('               - STT \n')
    end
end

% - Reflection Index (calcxulated using p1)
% pw_inds.RI_p1 = pw_inds.dia_amp./pw_inds.pulse_amp_p1;

% - Reflection Index (calculated using p2)
% pw_inds.RI_p2 = pw_inds.dia_amp./pw_inds.pulse_amp_p2;

% - Ratio of amplitudes of p2 and p1
% pw_inds.ratio_p2_p1 = pw_inds.pulse_amp_p2./pw_inds.pulse_amp_p1;

%% Areas and widths

if isfield(fid_pts, 'dic')
    do_area = 1;
else
    do_area = 0;
    if verbose_flag
        fprintf('               - Areas \n')
    end
end

store_width25 = nan(size(fid_pts.s.ind));
store_width50 = nan(size(fid_pts.s.ind));
store_width75 = nan(size(fid_pts.s.ind));


%params to reduce overhead
s = fid_pts.s;
ts = PPG.ts;
f1 = fid_pts.f1;
f2 = fid_pts.f2;
if do_area
    A1 = nan(size(fid_pts.s.ind));
    A2 = nan(size(fid_pts.s.ind));
    dic= fid_pts.dic;
else
    %preallocate
    A1 = [];
    A2 = [];
    dic= [];
end
fs = PPG.fs;
onsets = PPG.onsets;
t = PPG.t;


store_sVRI = nan(num_beats,1);

find_crossing = @(v, t) find((v(:)-t).*circshift((v(:)-t), [-1 0]) <= 0);

%%% Reduce the load of broadcast variables
f1_ind = f1.ind;
dic_ind = dic.ind;
f2_ind = f2.ind;
f1_amp = f1.amp;
s_amp = s.amp;



% parfor beat_no = 1 : length(PPG.onsets)-1
for beat_no = 1 : length(PPG.onsets)-1
    %     curr = [];
    %get current pulse
    if PPG.sqi_beat(beat_no) < 0.8
        continue        
    end
    
    curr_els = onsets(beat_no):onsets(beat_no+1);
    curr_ts = ts_filt(curr_els);
    correction_line = linspace(curr_ts(1), curr_ts(end), length(curr_ts));
    curr_ts = curr_ts - correction_line' + curr_ts(1);
    
    curr_t = t(curr_els);
    
    
    % skip if there isn't sufficient information for this pulse wave
    if isnan(s.ind(beat_no)) || isnan(f1_ind(beat_no)) || isnan(dic_ind(beat_no))
        continue
    end
    if do_area
        % find baseline of pulse wave (joining initial pulse onset to final pulse onset)
        % TO DO -- discuss do we use filtered or not filtered signal
        
        curr_dic_ind = dic_ind(beat_no) - f1_ind(beat_no) +1 ; 
        
        % - systolic area
        %calculate area -- sum up
%         A1(beat_no) = sum(curr_ts(1:curr_dic_ind))/( fs*pulse_amp(beat_no));
        A1(beat_no) = trapz(curr_ts(1:curr_dic_ind));
        
        mean_sys = mean(curr_ts(1:curr_dic_ind));
        
        % - diastolic area
        A2(beat_no) = trapz(curr_ts(curr_dic_ind:end));
        
        mean_dias = mean(curr_ts(curr_dic_ind:end));
    end
    
    store_sVRI(beat_no) = mean_dias/mean_sys;
    
    %%%% Widths

    
    %25
    cutoff_25 = 0.25*(s.amp(beat_no)-f1_amp(beat_no));
    cutoff_25 = cutoff_25 + f1_amp(beat_no);
    crossing_25 = find_crossing(curr_ts, cutoff_25);
    if length(crossing_25) < 2
        continue
    end
    store_width25(beat_no) = curr_t(crossing_25(end)) - curr_t(crossing_25(1));
    
    %50
    cutoff_50 = 0.5*(s.amp(beat_no)-f1_amp(beat_no));
    cutoff_50 = cutoff_50 + f1_amp(beat_no);
    crossing_50 = find_crossing(curr_ts, cutoff_50);
    if length(crossing_50) < 2
        continue
    end
    store_width50(beat_no) = curr_t(crossing_50(end)) - curr_t(crossing_50(1));
    
    
    %75
    cutoff_75 = 0.75*(s.amp(beat_no)-f1_amp(beat_no));
    cutoff_75 = cutoff_75 + f1_amp(beat_no);
    crossing_75 = find_crossing(curr_ts, cutoff_75);
    if length(crossing_75) < 2
        continue
    end
    store_width75(beat_no) = curr_t(crossing_75(end)) - curr_t(crossing_75(1));
    
    
end
if do_area
    pw_inds.A1 = A1;
    pw_inds.A2 = A2;
    
    % - Ratio of diastolic to systolic area (called Inflection point area)
    pw_inds.IPA = pw_inds.A2 ./ pw_inds.A1;
end
clear beat_no rel_beats


pw_inds.width25 = store_width25;
pw_inds.width50 = store_width50;
% pw_inds.width75 = store_width75; IGNORED as is affected by changing
% morphology at the top of the PPG 


%sVRI - Stress induced Vascular response index -- ﻿Lyu et al.
pw_inds.sVRI  = store_sVRI;

%% Frequency features

% - Normalised Harmonic Area (from Wang et al) - "Noninvasive cardiac output estimation using a novel PPG index"
% - Skewness and Kurtosis from Slapnicar et al
store_NHA = nan(num_beats,1);
store_skew  = nan(num_beats,1);
store_kurt  = nan(num_beats,1);

for beat_idx = 1:num_beats
    if PPG.sqi_beat(beat_idx) < sqi_threshold
        continue        
    end
    pulse = ts_filt(onsets(beat_idx):onsets(beat_idx+1));
    fft_opts.detrend = 1;
    pw = func.waveform.fft(pulse, PPG.fs, fft_opts);
    [~, loc] = findpeaks(pw);
    store_NHA(beat_idx) = 1-(sum(pw(loc(2:end)))/sum(pw(loc(1:end))));
    
    store_skew(beat_idx) = skewness(pulse);
    store_kurt(beat_idx) = kurtosis(pulse);
end
pw_inds.NHA = store_NHA;
pw_inds.skewness = store_skew;
pw_inds.kurtosis = store_kurt;

clear pw loc store_NHA

% - Inflection and Harmonic area ratio (IHAR) (from Wang et al) - "Noninvasive cardiac output estimation using a novel PPG index"
pw_inds.IHAR = pw_inds.NHA./pw_inds.IPA;



%% First Derivative

% - Maximum slope
% pw_inds.ms = fid_pts.amp.ms;

% - Maximum slope divided by the pulse amplitude
% pw_inds.ms_div_amp = fid_pts.amp.ms./pw_inds.pulse_amp;

% Get mean and variance of VPG in systolic and diastolic phase - Sun et al

pw_inds.sp_mean     = nan(num_beats,1);
pw_inds.sp_var      = nan(num_beats,1);
pw_inds.dp_mean     = nan(num_beats,1);
pw_inds.dp_var      = nan(num_beats,1);

for pulse_no = 1 : length(PPG.onsets)-1
    
    if PPG.sqi_beat(pulse_no) < sqi_threshold
        continue        
    end
    %get current pulse
    curr_els            = onsets(pulse_no):onsets(pulse_no+1);
    notch_loc           = PPG.fid_pts.dic.ind(pulse_no);
    sys_els = curr_els(curr_els < notch_loc);
    dia_els = curr_els(curr_els > notch_loc);
    
    sys_deriv        = PPG.derivs.first(sys_els);
    dia_deriv        = PPG.derivs.first(dia_els);
    
    pw_inds.sp_mean(pulse_no) = nanmean(sys_deriv);
    pw_inds.sp_var(pulse_no)  = nanvar(sys_deriv);
    pw_inds.dp_mean(pulse_no) = nanmean(dia_deriv);
    pw_inds.dp_var(pulse_no)  = nanvar(dia_deriv);
    
end

clear sys_deriv dia_deriv sys_els dia_els
%% Second Derivative

% - Amplitude of 'b' relative to 'a'
try
    pw_inds.b_div_a = fid_pts.b.amp./fid_pts.a.amp;
catch
    if verbose_flag
        fprintf('               - b/a \n')
    end
end

% - Amplitude of 'c' relative to 'a'
try
    pw_inds.c_div_a = fid_pts.c.amp./fid_pts.a.amp;
catch
    if verbose_flag
        fprintf('               - c/a \n')
    end
end

% - Amplitude of 'd' relative to 'a'
try
    pw_inds.d_div_a = fid_pts.d.amp./fid_pts.a.amp;
catch
    if verbose_flag
        fprintf('               - d/a \n')
    end
end

% - Amplitude of 'e' relative to 'a'
try
    pw_inds.e_div_a = fid_pts.e.amp./fid_pts.a.amp;
catch
    if verbose_flag
        fprintf('               - e/a \n')
    end
end

% % - Amplitude of 'a' relative to pulse amplitude
% try
%     pw_inds.a_div_amp = fid_pts.a.amp./pulse_amp;
% catch
%     if verbose_flag
%         fprintf('               - a/amp \n')
%     end
% end
% 
% % - Amplitude of 'b' relative to pulse amplitude
% try
%     pw_inds.b_div_amp = fid_pts.b.amp./pulse_amp;
% catch
%     if verbose_flag
%         fprintf('               - b/a \n')
%     end
% end
% 
% % - Amplitude of 'c' relative to pulse amplitude
% try
%     pw_inds.c_div_amp = fid_pts.c.amp./pulse_amp;
% catch
%     if verbose_flag
%         fprintf('               - b/amp \n')
%     end
% end
% 
% % - Amplitude of 'd' relative to pulse amplitude
% try
%     pw_inds.d_div_amp = fid_pts.d.amp./pulse_amp;
% catch
%     if verbose_flag
%         fprintf('               - b/amp \n')
%     end
% end
% 
% % - Amplitude of 'e' relative to pulse amplitude
% try
%     pw_inds.e_div_amp = fid_pts.e.amp./pulse_amp;
% catch
%     if verbose_flag
%         fprintf('               - b/amp \n')
%     end
% end

% - Ageing index: original
try
    pw_inds.AGI = pw_inds.b_div_a - pw_inds.c_div_a - pw_inds.d_div_a - pw_inds.e_div_a;
catch
    if verbose_flag
        fprintf('               - Ageing index \n')
    end
end

% - Ageing index: informal
% pw_inds.AGI_inf = pw_inds.b_div_a - pw_inds.e_div_a;

% - Ageing index: modified
% pw_inds.AGI_mod = pw_inds.b_div_a - pw_inds.c_div_a - pw_inds.d_div_a;

% Calculate SIs from second derivative slopes
try
    pt1.t = fid_pts.b.t;
    pt1.v = fid_pts.b.amp;
    pt2.t = fid_pts.c.t;
    pt2.v = fid_pts.c.amp;
    pw_inds.slope_b_c = ((pt2.v - pt1.v)./fid_pts.a.amp)./(pt2.t-pt1.t);
catch
    if verbose_flag
        fprintf('               - Slope b c \n')
    end
end

try
    pt2.t = fid_pts.d.t;
    pt2.v = fid_pts.d.amp;
    pw_inds.slope_b_d = ((pt2.v - pt1.v)./fid_pts.a.amp)./(pt2.t-pt1.t);
catch
    if verbose_flag
        fprintf('               - Slope b d \n')
    end
end


% - PPG AI - Pilt et al
try
    pw_inds.PPG_AI = nan(num_beats,1);
    notnan = ~isnan([fid_pts.d.ind + fid_pts.b.ind]);
    pw_inds.PPG_AI(notnan) = ts_filt(fid_pts.d.ind(notnan))./ts_filt(fid_pts.b.ind(notnan));
catch
    if verbose_flag
        fprintf('               - PPG AI \n')
    end
end

%% Pressure index
% - Pressure Index - Shin et al
try
    pw_inds.PI = (fid_pts.dic.t - fid_pts.s.t)./(fid_pts.dic.t - fid_pts.W.t) * ht;
catch
    if verbose_flag
        fprintf('               - Pressure index \n')
    end
end

%% Indices calculated from multiple derivatives

% - Ratio of diastolic to systolic area (called Inflection point area) plus d-peak
% pw_inds.IPAD = pw_inds.IPA + pw_inds.d_div_a;

% - Stiffness constant
% pw_inds.k = fid_pts.amp.s ./ ((fid_pts.amp.s - fid_pts.amp.ms ) ./ pw_inds.pulse_amp );

%% Gaussian decomposition
if do_gauss
    pw_inds.gauss_AI = nan(num_beats,1);
    pw_inds.gauss_RI = nan(num_beats,1);
    pw_inds.gauss_sys_dias = nan(num_beats,1);
    pw_inds.gauss_RTT_Rubins = nan(num_beats,1);
    pw_inds.gauss_AIx_Rubins = nan(num_beats,1);
    pw_inds.gauss_RI_Rubins = nan(num_beats,1);
    pw_inds.gauss_LVET = nan(num_beats,1);
    
    gaussian = @(b,x) b(1) * exp(-(x - b(2)).^2/b(3));
    
    for pulse_no = 1 : num_beats
        if PPG.sqi_beat(pulse_no) < sqi_threshold
            continue        
        end
        curr = [];
        %get current pulse
        curr_els = onsets(pulse_no):onsets(pulse_no+1);
        curr.t = PPG.t(curr_els);
        curr.t = curr.t - curr.t(1);
        %Get g1, g2, g3, g4
        for idx = 1:4
            eval(['g',num2str(idx),' = gaussian([fid_pts.g',num2str(idx),'.amp(pulse_no), fid_pts.g',num2str(idx),'.mu(pulse_no), fid_pts.g',num2str(idx),'.sigma(pulse_no)], curr.t);'])
        end
        systolic = g1 + g2;
        dias = g3+g4;
        
        
        pw_inds.gauss_AI(pulse_no) = max(systolic) - fid_pts.g3.amp(pulse_no);
        pw_inds.gauss_RI(pulse_no) = (sum(systolic) - sum(g3))./PPG.fs;
        pw_inds.gauss_sys_dias(pulse_no) = sum(systolic)/sum(dias);
        
        % Get LVET
        ds_dt = diff(systolic);
        ds2_dt2 = diff(ds_dt);
        ds3_dt3 = diff(ds2_dt2);
        
        
        [~, loc_peaks] = findpeaks(ds3_dt3);
        if length(loc_peaks) > 2
            pw_inds.gauss_LVET(pulse_no) = (loc_peaks(3) - loc_peaks(1))/fs;
        end
        
    end
    
    %Rubins et al
    pw_inds.gauss_RTT_Rubins = fid_pts.g3.mu - fid_pts.g1.mu; % Transit time of reflected wave
    pw_inds.gauss_AIx_Rubins = (fid_pts.g1.amp - fid_pts.g2.amp)./fid_pts.g1.amp; % Augmentation index
    pw_inds.gauss_RI_Rubins = fid_pts.g3.amp./fid_pts.g1.amp; % Reflection index
    
    %features I made from looking at videos
    pw_inds.gauss_amp4_amp1 = fid_pts.g4.amp./fid_pts.g1.amp;
    pw_inds.gauss_sigma4_amp1 = fid_pts.g4.sigma./fid_pts.g1.amp;

    %Add gaussian parameters as fiducial points
    pw_inds.g1_amp = fid_pts.g1.amp;
    pw_inds.g1_mu = fid_pts.g1.mu;
    pw_inds.g1_sigma = fid_pts.g1.sigma;
    pw_inds.g2_amp = fid_pts.g2.amp;
    pw_inds.g2_mu = fid_pts.g2.mu;
    pw_inds.g2_sigma = fid_pts.g2.sigma;
    pw_inds.g3_amp = fid_pts.g3.amp;
    pw_inds.g3_mu = fid_pts.g3.mu;
    pw_inds.g3_sigma = fid_pts.g3.sigma;
    pw_inds.g4_amp = fid_pts.g4.amp;
    pw_inds.g4_mu = fid_pts.g4.mu;
    pw_inds.g4_sigma = fid_pts.g4.sigma;

    
end


%% Set all when sqi < 0.8 to nan
set_sqi_func = @(x) set_sqi_pw(x, PPG.sqi_beat, sqi_threshold);
pw_inds = structfun(set_sqi_func, pw_inds, 'UniformOutput', false);
PPG.pw_inds = pw_inds;
%%
if plot_flag
    
    %if all data is nan then we cannot plot
    if sum(isnan(fid_pts.a.ind)) == length(fid_pts.a.ind)
        return
    end
    
    
    %plot a histogram of all pulse wave indicies
    
    
    pt_names_plt = fieldnames(pw_inds);
    
    colours = constants_def('Colours');
    
    %histogram of all points
    %First work ot number of rows and columns
    num_rows = floor(sqrt(numel(pt_names_plt)));
    num_columns = ceil(sqrt(numel(pt_names_plt)));
    while num_rows * num_columns < numel(pt_names_plt)
        if num_rows < num_columns
            num_rows = num_rows+1;
        else
            num_columns = num_columns+1;
        end
        
    end
    figure
    
    for pt_name_no = 1 : numel(pt_names_plt)
        subplot(num_rows, num_columns, pt_name_no)
        eval(['histogram(pw_inds.' pt_names_plt{pt_name_no},', ''FaceColor'', colours.blue)'])
        title_name =pt_names_plt{pt_name_no};
        title_name = strrep(title_name,'_', ' ');
        title(title_name)
    end
    
    %     func.plot.tightfig();
    %     set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 0.6);
    
end


end
%% Local funcs

function pw_out = set_sqi_pw(pw_in, sqi, sqi_threshold)
   pw_out = pw_in;
   pw_out(sqi < sqi_threshold) = nan;
end





