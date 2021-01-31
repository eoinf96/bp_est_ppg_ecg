% get_ppg_indices - This function computes the inicies from the PPG that
% are commonly used in BP estimation.
%Code adapted from P.Charlton: https://github.com/peterhcharlton/RRest

% Inputs : PPG -- This struct must contain all relevent fiducial points
% Outputs : pw_inds -- a struct of pulse wave indices computed for beat located
% in the PPG including the median value computed from beats of high signal
% quality
function [pw_inds] = get_norm_ppg_indices(PPG,  norm_ts, norm_derivs,  plot_flag,ht, sqi_threshold)
narginchk(1, inf);
PPG.norm_ts = norm_ts;
PPG.norm_derivs = norm_derivs;
if nargin < 4
    plot_flag = false;
end
if nargin < 5
    configs = constants_def('MOLLIE');
    volunteer_idx = PPG.record_name(2:3);
    volunteer_idx = str2double(volunteer_idx);
    ht = configs.demographics.height(volunteer_idx);
end
if nargin < 6
    sqi_threshold = 0.8;
end

if ~isfield(PPG, 'fid_pts')
    PPG.fid_pts= func.single_ppg.get_ppg_fid_pts(PPG);
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

% parfor beat_no = 1 : num_beats
for beat_no = 1 : length(PPG.onsets)-1
    
    % skip if there isn't sufficient information for this pulse wave
    if isnan(s.ind(beat_no)) || isnan(f1.ind(beat_no)) || isnan(dic.ind(beat_no))
        continue
    end
    
    curr = [];
    
    curr.ts = PPG.norm_ts{beat_no};
    curr.fs = length(curr.ts);
    curr.t = [0:(length(curr.ts) -1)]' /curr.fs;
    
    if curr.fs ==0
       continue 
    end
    
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
    
    
    
    
    
    
    %% Liang
    
    
    VPG = PPG.norm_derivs.first{beat_no};
    if ~any(isnan([fid_pts.b.ind(beat_no),fid_pts.d.ind(beat_no) , fid_pts.c.ind(beat_no) ]),2)
        
        store_liang_1(beat_no) = sum(curr.ts(fid_pts.f1.ind(beat_no):fid_pts.s.ind(beat_no)).^2);
        
        store_liang_2(beat_no) = (curr.ts(fid_pts.b.ind(beat_no)) - curr.ts(fid_pts.d.ind(beat_no)))./(fid_pts.b.t(beat_no) - fid_pts.d.t(beat_no));

        store_liang_4(beat_no) = curr.ts(fid_pts.c.ind(beat_no))./fid_pts.s.amp(beat_no);

        store_liang_8(beat_no) = VPG(fid_pts.c.ind(beat_no))./fid_pts.W.amp(beat_no);

        store_liang_10(beat_no) = (fid_pts.s.amp(beat_no) - curr.ts(fid_pts.c.ind(beat_no)))./(fid_pts.s.t(beat_no) - fid_pts.c.t(beat_no)); 
    end
    
       
    
    
    %% Frequency features
    
% - Normalised Harmonic Area (from Wang et al) - "Noninvasive cardiac output estimation using a novel PPG index"
% - Skewness and Kurtosis from Slapnicar et al
    curr.fs = length(curr.ts);
    if curr.fs ~=0
        [~, pw, ~, ~] = pe.sigproc.fft.fft(curr.ts, curr.fs, 1);
        [~, loc] = findpeaks(pw);

        store_NHA(beat_no) = 1-(sum(pw(loc(2:end)))/sum(pw(loc(1:end))));

        store_skew(beat_no) = skewness(curr.ts);
        store_kurt(beat_no) = kurtosis(curr.ts);


        
    end
     
    %% First derivative
        
    %get current pulse
    curr_els            = 1:length(PPG.norm_ts{beat_no});
    
    
    notch_loc           = fid_pts.dic.ind(beat_no);
    sys_els = curr_els(curr_els < notch_loc);
    dia_els = curr_els(curr_els > notch_loc);
    
    VPG = PPG.norm_derivs.first{beat_no};
    
    sys_deriv        = VPG(sys_els);
    dia_deriv        = VPG(dia_els);
    
    store_sp_mean(beat_no) = nanmean(sys_deriv);
    store_sp_var(beat_no)  = nanvar(sys_deriv);
    store_dp_mean(beat_no) = nanmean(dia_deriv);
    store_dp_var(beat_no)  = nanvar(dia_deriv);
    
    
    
    %% second derivative
    % - PPG AI - Pilt et al
     
   if ~or(isnan(fid_pts.d.ind(beat_no)), isnan(fid_pts.b.ind(beat_no)))
       store_PPG_AI(beat_no) = curr.ts(fid_pts.d.ind(beat_no))/curr.ts(fid_pts.b.ind(beat_no));
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
% 
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

% - Features from Liang, Elendi et al -- the top 10 features that correlate
% with BP
liang_numbers = 1:10;
for idx = liang_numbers
    eval(['pw_inds.liang_',num2str(idx),' = nan(num_beats,1);'])
end

pw_inds.liang_1 = store_liang_1;

pw_inds.liang_2 = store_liang_2;


% 3 - Time between S and c wave
pw_inds.liang_3 = fid_pts.s.t - fid_pts.c.t;

% 5 - Time between S and d
pw_inds.liang_5 = fid_pts.s.t - fid_pts.d.t;

pw_inds.liang_4 = store_liang_4;


% 6 - Ageing index
pw_inds.liang_6 = (fid_pts.b.amp - fid_pts.c.amp -fid_pts.d.amp)./(fid_pts.a.amp);

% 7 - Amplitude of d
pw_inds.liang_7 = fid_pts.d.amp;

pw_inds.liang_8 = store_liang_8;

% 9 - Ratio of d to a
pw_inds.liang_9 = fid_pts.d.amp./fid_pts.a.amp;


pw_inds.liang_10 = store_liang_10;


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
            pw_inds.gauss_LVET(pulse_no) = (loc_peaks(3) - loc_peaks(1))/( fs*pulse_amp(pulse_no));
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

