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
function [pw_inds] = get_norm_ppg_indices(PPG, plot_flag,database,do_filter, sqi_threshold)
%get_ppg_indices - This function computes the indicies from the PPG that
% are commonly used in BP estimation.
%Code adapted from P.Charlton: https://github.com/peterhcharlton/RRest

% Inputs : PPG -- This struct must contain all relevent fiducial points
% Outputs : pw_inds -- a struct of pulse wave indices computed for beat located
% in the PPG including the median value computed from beats of high signal
% quality
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
        %update this when the study is properly run
        ht = configs.demographics.height;
    end
end
if nargin < 4
    do_filter = 1;
end
if nargin < 5
    sqi_threshold = 0.8;
end

if ~isfield(PPG, 'norm_fid_pts')
    PPG.norm_fid_pts= func.pulsew.get_ppg_fid_pts(PPG);
end
fid_pts = PPG.norm_fid_pts;

do_gauss = 0;
if isfield(PPG.norm_fid_pts, 'g1')
    do_gauss = 1;
end
verbose_flag = 0;
fs = PPG.fs;

onsets = PPG.onsets;

num_beats = length(PPG.peaks);

%% Filter Signal
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
s_g_filter_len = 7;
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

%% Amplitudes

try
    pw_inds.dic_amp = fid_pts.dic.amp;
catch
    if verbose_flag
        fprintf('               - Dicrotic notch amplitude \n')
    end
end



% - Reflection Index 
try
    pw_inds.RI = fid_pts.dia.amp;
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


%% Do looped features

pulse_amp= nan(size(fid_pts.s.ind));

% Areas and widths
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
find_crossing = @(v, t) find((v(:)-t).*circshift((v(:)-t), [-1 0]) <= 0);

store_sVRI = nan(num_beats,1);

for idx = 1:10
    eval(['store_liang_',num2str(idx),' = nan(num_beats,1);'])
end


store_NHA = nan(num_beats,1);
store_skew  = nan(num_beats,1);
store_kurt  = nan(num_beats,1);

% Get mean and variance of VPG in systolic and diastolic phase - Sun et al

store_sp_mean     = nan(num_beats,1);
store_sp_var      = nan(num_beats,1);
store_dp_mean     = nan(num_beats,1);
store_dp_var      = nan(num_beats,1);


% Second derivative feature
store_PPG_AI = nan(num_beats,1);


% Gaussian features -- have a rethink here
if do_gauss
    pw_inds.gauss_AI = nan(num_beats,1);
    pw_inds.gauss_RI = nan(num_beats,1);
    pw_inds.gauss_sys_dias = nan(num_beats,1);
    pw_inds.gauss_RTT_Rubins = nan(num_beats,1);
    pw_inds.gauss_AIx_Rubins = nan(num_beats,1);
    pw_inds.gauss_RI_Rubins = nan(num_beats,1);
    pw_inds.gauss_LVET = nan(num_beats,1);
    gaussian = @(b,x) b(1) * exp(-(x - b(2)).^2/b(3));
end

%



% parfor beat_no = 1 : num_beats
for beat_no = 1 : length(PPG.onsets)-1
    
    % skip if there isn't sufficient information for this pulse wave
    if isnan(s.ind(beat_no)) || isnan(f1.ind(beat_no)) || isnan(dic.ind(beat_no))
        continue
    end
    
    if PPG.sqi_beat(beat_no) <sqi_threshold
        continue        
    end
    %% Get curr
    curr_els = onsets(beat_no):onsets(beat_no+1);
    
    curr = [];
    curr.t = PPG.t(curr_els);
    curr.t  = curr.t - curr.t(1);
    curr.fs = length(curr.t);
    curr.t = curr.t / curr.t(end);
    
    dt = 1/curr.fs;
    
    curr.ts = ts_filt(curr_els);
    curr.ts = curr.ts - min(curr.ts);
    curr.ts = curr.ts/ max(curr.ts);
    
    
    % Correct for low frequency baseline drift in a single beat
    correction_line = linspace(curr.ts(1), curr.ts(end), length(curr.ts));
    curr.ts = curr.ts - correction_line' + curr.ts(1);
    
    clear correction_line
    
    %
    curr.derivs.first   = func.waveform.savitzky_golay_deriv(curr.ts, 1, s_g_filter_len)./dt;     
    if curr.fs ==0
       continue 
    end
    
    %%
    
    pulse_amp(beat_no) = max(curr.ts);
    
    if do_area
        % find baseline of pulse wave (joining initial pulse onset to final pulse onset)

        baseline = linspace(curr.ts(f1.ind(beat_no)),...      %start
            curr.ts(f2.ind(beat_no)),...                      %end
            f2.ind(beat_no) - f1.ind(beat_no)+1 ...  %number of points
            );
        baseline = baseline(:);
        
        
        % - systolic area
        rel_pts = f1.ind(beat_no) : dic.ind(beat_no);
        
        if ~isnan(rel_pts)
            baseline_pts = rel_pts - f1.ind(beat_no) + 1;
            %calculate area -- sum up
            A1(beat_no) = sum(curr.ts(rel_pts) - baseline(baseline_pts))/( curr.fs);
            
            mean_sys = mean(curr.ts(rel_pts));

            % - diastolic area
            rel_pts = dic.ind(beat_no) : f2.ind(beat_no);
            baseline_pts = rel_pts - f1.ind(beat_no) + 1;
            A2(beat_no) = sum(curr.ts(rel_pts) - baseline(baseline_pts))/(curr.fs);
            
            mean_dias = mean(curr.ts(rel_pts));
        end
    end
    
    
    store_sVRI(beat_no) = mean_dias/mean_sys;
    
    
    %25
    cutoff_25 = 0.25*(s.amp(beat_no)-f1.amp(beat_no));
    cutoff_25 = cutoff_25 + f1.amp(beat_no);
    crossing_25 = find_crossing(curr.ts, cutoff_25);
    if ~length(crossing_25) < 2
        store_width25(beat_no) = curr.t(crossing_25(end)) - curr.t(crossing_25(1));
    end
    
    
    %50
    cutoff_50 = 0.5*(s.amp(beat_no)-f1.amp(beat_no));
    cutoff_50 = cutoff_50 + f1.amp(beat_no);
    crossing_50 = find_crossing(curr.ts, cutoff_50);
    if ~length(crossing_50) < 2
        store_width50(beat_no) = curr.t(crossing_50(end)) - curr.t(crossing_50(1));
    end
    
    
    
    %75
    cutoff_75 = 0.75*(s.amp(beat_no)-f1.amp(beat_no));
    cutoff_75 = cutoff_75 + f1.amp(beat_no);
    crossing_75 = find_crossing(curr.ts, cutoff_75);
    if ~length(crossing_75) < 2
        store_width75(beat_no) = curr.t(crossing_75(end)) - curr.t(crossing_75(1));
    end
    
    
    
    
    %% Frequency features
    
% - Normalised Harmonic Area (from Wang et al) - "Noninvasive cardiac output estimation using a novel PPG index"
% - Skewness and Kurtosis from Slapnicar et al
    if curr.fs ~=0
        [~, pw, ~, ~] = pe.sigproc.fft.fft(curr.ts, curr.fs, 1);
        [~, loc] = findpeaks(pw);

        store_NHA(beat_no) = 1-(sum(pw(loc(2:end)))/sum(pw(loc(1:end))));

        store_skew(beat_no) = skewness(curr.ts);
        store_kurt(beat_no) = kurtosis(curr.ts);
    end
     
    %% First derivative
            
    notch_loc           = fid_pts.dic.ind(beat_no);    
    VPG = curr.derivs.first;
    
    sys_deriv        = VPG(1:notch_loc);
    dia_deriv        = VPG(notch_loc:end);
    
    store_sp_mean(beat_no) = mean(sys_deriv, 'omitnan');
    store_sp_var(beat_no)  = var(sys_deriv, 'omitnan');
    store_dp_mean(beat_no) = mean(dia_deriv, 'omitnan');
    store_dp_var(beat_no)  = var(dia_deriv, 'omitnan');
    
    
    
    %% second derivative
    % - PPG AI - Pilt et al
     
   if ~or(isnan(fid_pts.d.ind(beat_no)), isnan(fid_pts.b.ind(beat_no)))
       store_PPG_AI(beat_no) = curr.ts(fid_pts.d.ind(beat_no))/curr.ts(fid_pts.b.ind(beat_no));
   end
   
   %% Do Gauss features here
   
   if do_gauss 
        %Get g1, g2, g3, g4
        for idx = 1:4
            eval(['g',num2str(idx),' = gaussian([fid_pts.g',num2str(idx),'.amp(beat_no), fid_pts.g',num2str(idx),'.mu(beat_no), fid_pts.g',num2str(idx),'.sigma(beat_no)], curr.t);'])
        end
        systolic = g1 + g2;
        dias = g3+g4;
        
        pw_inds.gauss_AI(beat_no) = max(systolic) - fid_pts.g3.amp(beat_no);
        pw_inds.gauss_RI(beat_no) = (sum(systolic) - sum(g3))./PPG.fs;
        pw_inds.gauss_sys_dias(beat_no) = sum(systolic)/sum(dias);
        
        % Get LVET
        ds_dt = diff(systolic);
        ds2_dt2 = diff(ds_dt);
        ds3_dt3 = diff(ds2_dt2);
        
        
        [~, loc_peaks] = findpeaks(ds3_dt3);
        if length(loc_peaks) > 2
            pw_inds.gauss_LVET(beat_no) = (loc_peaks(3) - loc_peaks(1))/( fs*pulse_amp(beat_no));
        end
   
   end
    
end
curr = []; pw = [];
if do_area
    pw_inds.A1 = A1;
    pw_inds.A2 = A2;
    
    % - Ratio of diastolic to systolic area (called Inflection point area)
    pw_inds.IPA = pw_inds.A2 ./ pw_inds.A1;
end
clear beat_no rel_beats 

%sVRI - Stress induced Vascular response index -- ﻿Lyu et al.
pw_inds.sVRI  = store_sVRI;


pw_inds.width25 = store_width25;
pw_inds.width50 = store_width50;
% pw_inds.width75 = store_width75; IGNORED as is affected by changing
% morphology at the top of the PPG 


pw_inds.NHA = store_NHA;
pw_inds.skewness = store_skew;
pw_inds.kurtosis = store_kurt;

clear pw loc store_NHA

% - Inflection and Harmonic area ratio (IHAR) (from Wang et al) - "Noninvasive cardiac output estimation using a novel PPG index"
pw_inds.IHAR = pw_inds.NHA./pw_inds.IPA;



pw_inds.sp_mean = store_sp_mean;
pw_inds.sp_var = store_sp_var;
pw_inds.dp_mean = store_dp_mean;
pw_inds.dp_var = store_dp_var;
clear sys_deriv dia_deriv sys_els dia_els




pw_inds.PPG_AI = store_PPG_AI;

%% Remaining Gauss features
if do_gauss
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

% % - Ageing index: original
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


