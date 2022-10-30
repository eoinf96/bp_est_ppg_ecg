function [pw_inds] = get_ppg_indices(PPG,do_normalise, config,  plot_flag)
%  This function computes the inicies from the PPG that are commonly used in BP estimation.
% Code adapted from P.Charlton:
% https://github.com/peterhcharlton/pulse-analyse
%
% INPUT: PPG -- This struct must contain all relevent fiducial points
%       do_normalise: whether to produce pw inds on normalised beats 
%       configs: configs struct, see below for details
%       plot_flag: flag whether to plot a summary
% 
% OUTPUT: pw_inds -- a struct of pulse wave indices computed for beat located in the PPG 
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
% - Padilla, J. M. et al. Assessment of relationships between blood pressure, pulse wave velocity and digital volume pulse. In Computers in Cardiology, vol. 33, 893–896 (IEEE, 2006).
% - Mok Ahn, J. New aging index using signal features of both photoplethysmograms and acceleration plethysmograms. Healthc. Informatics Res. 23, 53–59, DOI: 10.4258/hir.2017.23.1.53 (2017).
% - Elgendi, M. On the Analysis of Fingertip Photoplethysmogram Signals. Curr. Cardiol. Rev. 8, 14–25, DOI: 10.2174/157340312801215782 (2012).
% - Addison, P. S. Slope transit time (STT): A pulse transit time proxy requiring only a single signal fiducial point. IEEE Transactions on Biomed. Eng. 63, 2441–2444 (2016).
% - Lyu, Y. et al. Measuring photoplethysmogram-based stress-induced vascular response index to assess cognitive load and stress. In Proceedings of the 33rd annual ACM conference on human factors in computing systems, 857–866 (2015).
% - Wang, L., Pickwell-MacPherson, E., Liang, Y. P. & Zhang, Y. T. Noninvasive cardiac output estimation using a novel photoplethysmogram index. In Proceedings of the 31st Annual International Conference of the IEEE Engineering in Medicine and Biology Society: Engineering the Future of Biomedicine, EMBC 2009, 1746–1749, DOI: 10.1109/IEMBS.2009.5333091 (IEEE, 2009).
% - Rubins, U., Grabovskis, A., Grube, J. & Kukulis, I. Photoplethysmography analysis of artery properties in patients with cardiovascular diseases. In IFMBE Proceedings, vol. 20 IFMBE, 319–322, DOI: 10.1007/978-3-540-69367-3-85 (Springer, 2008)
narginchk(1, inf);
if nargin < 2 || isempty(do_normalise)
    do_normalise = 1; % Return PPG indicies of a normalised PPG beat.
end
if nargin < 3 || isempty(config)
    config = struct();
end
if nargin < 4
    plot_flag  = false;
end
default_config.do_normalise = do_normalise; % Whether to also return the normalised pulse wave inds
default_config.ht = nan; % The height of the individual
default_config.do_filter = 1; % Flags whether to filter the PPG signal
default_config.ht = nan; % Height of the individual
default_config.verbose_flag = 0;
config = func.aux_functions.update_with_default_opts(config, default_config);
%% Get fid_pts
if ~isfield(PPG, 'fid_pts')
    if config.do_normalise ~= 1
        [fid_pts, PPG.derivs] = func.pulsew.get_ppg_fid_pts(PPG);
    else
        [~, fid_pts, ~] = func.pulsew.get_ppg_fid_pts(PPG);
    end
else
    if config.do_normalise ~= 1
        fid_pts = PPG.fid_pts;
    else
        fid_pts = PPG.norm_fid_pts;
    end
    
    if ~isfield(PPG, 'derivs')
        PPG.derivs = func.pulsew.get_derivs(PPG.ts, PPG.fs);
    end
end

if isempty(fid_pts)
    pw_inds = [];
    return
end
num_beats = length(PPG.peaks);
%% Filter Signal
if config.do_filter    
    % Butterworth IIR bandpass filter
    [b,a] = butter(8,10/(PPG.fs/2), 'low');
    ts_filt = filtfilt(b,a,PPG.ts);
else
    ts_filt = PPG.ts;
end


%% Timings
if config.verbose_flag
    fprintf('Cannot get ...\n')
end

% - Time between systolic and diastolic peaks (secs)
try
    pw_inds.delta_t = fid_pts.dia.t-fid_pts.s.t;
catch
    if config.verbose_flag
        fprintf('               - Delta T \n')
    end
end

% - Crest time (secs)
try
    pw_inds.CT = fid_pts.s.t-fid_pts.f1.t;
catch
    if config.verbose_flag
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
    if config.verbose_flag
        fprintf('               - Duration of systole \n')
    end
end

% - Duration of diastole
try
    pw_inds.t_dia = fid_pts.f2.t-fid_pts.dic.t;
catch
    if config.verbose_flag
        fprintf('               - Duration of diastole \n')
    end
end

% - Ratio of systolic to diastolic durations
try
    pw_inds.t_ratio = pw_inds.t_sys./pw_inds.t_dia;
catch
    if config.verbose_flag
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
% try
%     pulse_amp = fid_pts.s.amp-fid_pts.f1.amp;
% catch
%     if config.verbose_flag
%         fprintf('               - Pulse Amplitude \n')
%     end
% end

% - Augmentation Pressure
% try
%     pw_inds.AP = fid_pts.p2pk.amp-fid_pts.p1in.amp;
% catch
%     if config.verbose_flag
%         fprintf('               - Augmentation Pressure \n')
%     end
% end


% % - Agumentation Index Takazawa
% try
%     pw_inds.AI = PPG.fid_pts.dia.amp./PPG.fid_pts.s.amp;
% catch
%     if config.verbose_flag
%         fprintf('               - Augmentation Index \n')
%     end
% end

% - Dicrotic notch amplitude
try
    pw_inds.dic_amp = fid_pts.dic.amp;
catch
    if config.verbose_flag
        fprintf('               - Dicrotic notch amplitude \n')
    end
end

% - Reflection Index (calculated using systolic peak)
try
    pw_inds.RI = (fid_pts.dia.amp - fid_pts.f1.amp);
catch
    if config.verbose_flag
        fprintf('               - Reflection index \n')
    end
end

% - Slope Transit Time (from Addison et al)
try
    pw_inds.STT = (fid_pts.s.t - fid_pts.f1.t)./(fid_pts.s.amp - fid_pts.f1.amp);
catch
    if config.verbose_flag
        fprintf('               - STT \n')
    end
end

% - Reflection Index (calcxulated using p1)
% pw_inds.RI_p1 = pw_inds.dia_amp./pw_inds.pulse_amp_p1;

% - Reflection Index (calculated using p2)
% pw_inds.RI_p2 = pw_inds.dia_amp./pw_inds.pulse_amp_p2;

% - Ratio of amplitudes of p2 and p1
% pw_inds.ratio_p2_p1 = pw_inds.pulse_amp_p2./pw_inds.pulse_amp_p1;


%% Second Derivative


% - Amplitude of 'b' relative to 'a'
try
    pw_inds.b_div_a = fid_pts.b.amp./fid_pts.a.amp;
catch
    if config.verbose_flag
        fprintf('               - b/a \n')
    end
end

% - Amplitude of 'c' relative to 'a'
try
    pw_inds.c_div_a = fid_pts.c.amp./fid_pts.a.amp;
catch
    if config.verbose_flag
        fprintf('               - c/a \n')
    end
end

% - Amplitude of 'd' relative to 'a'
try
    pw_inds.d_div_a = fid_pts.d.amp./fid_pts.a.amp;
catch
    if config.verbose_flag
        fprintf('               - d/a \n')
    end
end

% - Amplitude of 'e' relative to 'a'
try
    pw_inds.e_div_a = fid_pts.e.amp./fid_pts.a.amp;
catch
    if config.verbose_flag
        fprintf('               - e/a \n')
    end
end

% % - Ageing index: original
try
    pw_inds.AGI = pw_inds.b_div_a - pw_inds.c_div_a - pw_inds.d_div_a - pw_inds.e_div_a;
catch
    if config.verbose_flag
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
    if config.verbose_flag
        fprintf('               - Slope b c \n')
    end
end

try
    pt2.t = fid_pts.d.t;
    pt2.v = fid_pts.d.amp;
    pw_inds.slope_b_d = ((pt2.v - pt1.v)./fid_pts.a.amp)./(pt2.t-pt1.t);
catch
    if config.verbose_flag
        fprintf('               - Slope b d \n')
    end
end



%% Pressure index
% - Pressure Index - Shin et al
try
    pw_inds.PI = (fid_pts.dic.t - fid_pts.s.t)./(fid_pts.dic.t - fid_pts.W.t) * config.ht;
catch
    if config.verbose_flag
        fprintf('               - Pressure index \n')
    end
end

%% Indices calculated from multiple derivatives

% - Ratio of diastolic to systolic area (called Inflection point area) plus d-peak
% pw_inds.IPAD = pw_inds.IPA + pw_inds.d_div_a;

% - Stiffness constant
% pw_inds.k = fid_pts.amp.s ./ ((fid_pts.amp.s - fid_pts.amp.ms ) ./ pw_inds.pulse_amp );

%% Do looped features

% parfor beat_no = 1 : length(PPG.onsets)-1
for beat_no = 1 : num_beats
    %     curr = [];
    %get current pulse
    if PPG.sqi_beat(beat_no) ==0 
        continue        
    end
    
    curr_els = PPG.onsets(beat_no):PPG.onsets(beat_no+1);
    curr.ts = ts_filt(curr_els);
    correction_line = linspace(curr.ts(1), curr.ts(end), length(curr.ts));
    curr.ts = curr.ts - correction_line' + curr.ts(1);
    
    curr.t = PPG.t(curr_els);curr.t = curr.t - curr.t(1);
    curr.fs = PPG.fs;
    
    curr.derivs.first = PPG.derivs.first(curr_els);
    
    if config.do_normalise
        curr.ts = curr.ts - min(curr.ts);
        curr.ts = curr.ts/ max(curr.ts);
        curr.fs = length(curr.t);
        curr.t = curr.t / curr.t(end);
        % Get derivatives
        curr.derivs = func.pulsew.get_derivs(curr.ts, curr.fs);
    end    
    pw_inds = local_get_pulse_inds(fid_pts, curr, pw_inds, beat_no);    
end

 % - Ratio of diastolic to systolic area (called Inflection point area)
pw_inds.IPA = pw_inds.A2 ./ pw_inds.A1;

% - Inflection and Harmonic area ratio (IHAR) (from Wang et al) - "Noninvasive cardiac output estimation using a novel PPG index"
pw_inds.IHAR = pw_inds.NHA./pw_inds.IPA;

%% Remaining Gauss features
if isfield(fid_pts, 'g1')
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

%% make all indices be a column vector
pw_names = fieldnames(pw_inds);
for p_idx = 1:length(pw_names)
    if all(isnan(pw_inds.(pw_names{p_idx})))
        pw_inds = rmfield(pw_inds, pw_names{p_idx});
    else
        pw_inds.(pw_names{p_idx}) = pw_inds.(pw_names{p_idx})(:);
    end
end

%% PLOT
if plot_flag
    %if all data is nan then we cannot plot
    if sum(isnan(fid_pts.a.ind)) == length(fid_pts.a.ind)
        return
    end
    %plot a histogram of all pulse wave indicies
    pt_names_plt = fieldnames(pw_inds);
    
    colours = func.aux_functions.define_colours;
    
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
    figure('position', [ 425         191        1496         778])
    
    for pt_name_no = 1 : numel(pt_names_plt)
        subplot(num_rows, num_columns, pt_name_no)
        histogram(pw_inds.(pt_names_plt{pt_name_no}), 'FaceColor', colours.grey)
        title(strrep(pt_names_plt{pt_name_no},'_', ' '))
        ylabel('#')
    end
    func.plot.tightfig();
end


end

%% Local funcs

function pw_inds = local_get_pulse_inds(fid_pts, curr, pw_inds, beat_no)
    % Local function to return all pw inds that require assessment on a beat by beat level
    find_crossing = @(v, t) find((v(:)-t).*circshift((v(:)-t), [-1 0]) <= 0);

    
    % Skip if there isn't sufficient information for this pulse wave
    if isnan(fid_pts.s.ind(beat_no)) || isnan(fid_pts.f1.ind(beat_no)) || isnan(fid_pts.dic.ind(beat_no))
        return
    end
    
    starting_ind = fid_pts.f1.ind(beat_no) - 1;
    %% AREAS
    % find baseline of pulse wave (joining initial pulse onset to final pulse onset)
    curr_dic_ind = fid_pts.dic.ind(beat_no) - starting_ind ; 

    % - systolic area
    %calculate area -- sum up
    %         A1(beat_no) = sum(curr.ts(1:curr_dic_ind))/( fs*pulse_amp(beat_no));
    pw_inds.A1(beat_no) = trapz(curr.ts(1:curr_dic_ind))/curr.fs;

    mean_sys = mean(curr.ts(1:curr_dic_ind));

    % - diastolic area
    pw_inds.A2(beat_no) = trapz(curr.ts(curr_dic_ind:end))/curr.fs;

    mean_dias = mean(curr.ts(curr_dic_ind:end));
        
    pw_inds.sVRI(beat_no) = mean_dias/mean_sys;
    
    %% Widths 
    %25
    cutoff_25 = 0.25*(fid_pts.s.amp(beat_no)-fid_pts.f1.amp(beat_no));
    cutoff_25 = cutoff_25 + fid_pts.f1.amp(beat_no);
    crossing_25 = find_crossing(curr.ts, cutoff_25);
    if length(crossing_25) >1
        pw_inds.width25(beat_no) = curr.t(crossing_25(end)) - curr.t(crossing_25(1));
    end
    
    
    %50
    cutoff_50 = 0.5*(fid_pts.s.amp(beat_no)-fid_pts.f1.amp(beat_no));
    cutoff_50 = cutoff_50 + fid_pts.f1.amp(beat_no);
    crossing_50 = find_crossing(curr.ts, cutoff_50);
    if length(crossing_50) >1
        pw_inds.width50(beat_no) = curr.t(crossing_50(end)) - curr.t(crossing_50(1));
    end
    
    
    
    %75
    cutoff_75 = 0.75*(fid_pts.s.amp(beat_no)-fid_pts.f1.amp(beat_no));
    cutoff_75 = cutoff_75 + fid_pts.f1.amp(beat_no);
    crossing_75 = find_crossing(curr.ts, cutoff_75);
    if length(crossing_75) >1
        pw_inds.width75(beat_no) = curr.t(crossing_75(end)) - curr.t(crossing_75(1));
    end
    
    
    %% Frequency features
    
% - Normalised Harmonic Area (from Wang et al) - "Noninvasive cardiac output estimation using a novel PPG index"
% - Skewness and Kurtosis from Slapnicar et al
    if curr.fs ~=0
        fft_opts.detrend = 1;
        pw = func.waveform.fft(curr.ts, curr.fs, fft_opts);
        [~, loc] = findpeaks(pw);

        pw_inds.NHA(beat_no) = 1-(sum(pw(loc(2:end)))/sum(pw(loc(1:end))));

        pw_inds.skewness(beat_no) = skewness(curr.ts);
        pw_inds.kurtosis(beat_no) = kurtosis(curr.ts);
    end
     
    %% First derivative
            
    notch_loc           = fid_pts.dic.ind(beat_no) - starting_ind;    
    VPG = curr.derivs.first;
    
    sys_deriv        = VPG(1:notch_loc);
    dia_deriv        = VPG(notch_loc:end);
    
    pw_inds.sp_mean(beat_no) = mean(sys_deriv, 'omitnan');
    pw_inds.sp_var(beat_no)  = var(sys_deriv, 'omitnan');
    pw_inds.dp_mean(beat_no) = mean(dia_deriv, 'omitnan');
    pw_inds.dp_var(beat_no)  = var(dia_deriv, 'omitnan');
   
    %% second derivative
    % - PPG AI - Pilt et al
     
   if ~or(isnan(fid_pts.d.ind(beat_no)), isnan(fid_pts.b.ind(beat_no)))
       pw_inds.PPG_AI(beat_no) = curr.ts(fid_pts.d.ind(beat_no) - starting_ind)/curr.ts(fid_pts.b.ind(beat_no) - starting_ind);
   end
   
   %% Do Gauss features here
   
   if isfield(fid_pts, 'g1')
       gaussian = @(b,x) b(1) * exp(-(x - b(2)).^2/b(3));
 
        %Get g1, g2, g3, g4
        for idx = 1:4
            eval(['g',num2str(idx),' = gaussian([fid_pts.g',num2str(idx),'.amp(beat_no), fid_pts.g',num2str(idx),'.mu(beat_no), fid_pts.g',num2str(idx),'.sigma(beat_no)], curr.t);'])
        end
        systolic = g1 + g2;
        dias = g3+g4;
        
        pw_inds.gauss_AI(beat_no) = max(systolic) - fid_pts.g3.amp(beat_no);
        pw_inds.gauss_RI(beat_no) = (sum(systolic) - sum(g3))./curr.fs;
        pw_inds.gauss_sys_dias(beat_no) = sum(systolic)/sum(dias);
        
        % Get LVET
        ds_dt = diff(systolic);
        ds2_dt2 = diff(ds_dt);
        ds3_dt3 = diff(ds2_dt2);
        
        
        [~, loc_peaks] = findpeaks(ds3_dt3);
        if length(loc_peaks) > 2
            pw_inds.gauss_LVET(beat_no) = (loc_peaks(3) - loc_peaks(1))/( curr.fs);
        end
   
   end
end

