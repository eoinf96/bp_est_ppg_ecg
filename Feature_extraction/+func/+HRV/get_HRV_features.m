function [HRV_stats] = get_HRV_features(RR_ts, RR_times, t_window_start, t_window_end, configs)
%  This function returns a feature vector formed from ECG features typically
%  used in BP estimation
%
% INPUT: RR_ts: RR intervals
%        RR_ts RR intervals times
%        t_window_start: vector of window start times within which ECG features are computed.
%        t_window_end: vector of window end times within which ECG features are computed.
%        configs: configs struct, see below for details
%
% OUTPUT: HRV_stats: Struct of HRV features
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
% Relevant literature:


narginchk(2, inf);
if isempty(RR_times) % Return a nothing
    HRV_stats.SDRR = nan(1,1);
    HRV_stats.pRR50 = nan(1,1); %Measure of parasympathetic activity
    HRV_stats.RMSSD = nan(1,1);
    
    %Frequency domain
    HRV_stats.VLF_power = nan(1,1);
    HRV_stats.nLF_power = nan(1,1);
    HRV_stats.nHF_power = nan(1,1);
    HRV_stats.T_power   = nan(1,1);
    HRV_stats.LF_div_HF = nan(1,1);
    
    % Non linear
    HRV_stats.poincare_SDD1 = nan(1,1); %thought to be equivalent to RMSSD - short term HRV (sympathetic?)%Vasoconstriction is influenced by the low frequency sympathetic nervous system
    HRV_stats.poincare_SDD2 = nan(1,1);
    HRV_stats.poincare_area = nan(1,1);
    HRV_stats.poincare_SDD1_SDD2 = nan(1,1);
    return
end

if nargin < 3 || isempty(t_window_start) || isempty(t_window_end)
    %Recommended by the HR variability task force Camm et al 1996
    window_length_time = 60;
    window_step_time = 60;
    t_window_start = 0:window_step_time:(RR_times(end)-window_step_time);
    t_window_end = t_window_start + window_length_time;
end
% Set default params
if nargin < 5 || isempty(configs)
   configs = struct(); 
end
default_configs.plot_flag = 0;
default_configs.signal_name = 'RR interval';
configs = func.aux_functions.update_with_default_opts(configs, default_configs);


if length(RR_ts) ~= length(RR_times)
    error('RR int and times have to be the same length')
end
%% Parameter definitions
%define frequency bands
VLFP_low = 0.0033; VLFP_high = 0.04;
LFP_low = 0.04;LFP_high = 0.16; % Sympathetic + parasympathetic (although debated)
HFP_low = 0.16;HFP_high = 0.40; % Parasympathetic activity and Respiration
%% Interpolate RR signal so that we can compute an accurate estimate of its frequency spectrum
fs_new=1;

RR.ts_beat = RR_ts;
RR.t_beat = RR_times;
RR.t = 0:(1/fs_new):RR.t_beat(end);
RR.ts = interp1( RR.t_beat, RR.ts_beat,  RR.t, 'PCHIP', 'extrap');
%% Define the window start and end times
if isempty(t_window_start)
    t_window_start = RR_times(1);
    t_window_end = RR_times(end);
end

window_start = floor(t_window_start *fs_new)+1;
window_end = floor(t_window_end *fs_new);
window_end(window_end > length(RR.t)) = length(RR.t);
num_windows = length(t_window_start);


%% Initialise features
%Time domain

HRV_stats.SDRR = nan(num_windows,1);
HRV_stats.pRR50 = nan(num_windows,1); %Measure of parasympathetic activity
HRV_stats.RMSSD = nan(num_windows,1);

%Frequency domain
HRV_stats.VLF_power = nan(num_windows,1);
HRV_stats.nLF_power = nan(num_windows,1);
HRV_stats.nHF_power = nan(num_windows,1);
HRV_stats.T_power   = nan(num_windows,1);
% HRV_stats.LF_div_HF = nan(num_windows,1);

% Non linear
HRV_stats.poincare_SDD1 = nan(num_windows,1); %thought to be equivalent to RMSSD - short term HRV (sympathetic?)
HRV_stats.poincare_SDD2 = nan(num_windows,1);
% poincare_area = nan(num_windows,1);

%Rotation matrix for use in poincare analysis
R = [cos(-pi/4) -sin(-pi/4);sin(-pi/4)  cos(-pi/4)] ;


for wind_idx = 1:num_windows
    hrv_window_beat = RR.ts_beat(and(RR.t_beat> t_window_start(wind_idx),RR.t_beat < t_window_end(wind_idx) ));
    hrv_window = RR.ts(window_start(wind_idx):window_end(wind_idx));
    if length(hrv_window) < 2
        continue
    end
    %% Time domain
    HRV_stats.SDRR(wind_idx) = std(hrv_window_beat);
    HRV_stats.pRR50(wind_idx) = sum(abs(diff(hrv_window_beat)) > 0.050 )/(length(hrv_window_beat)-1);
    HRV_stats.RMSSD(wind_idx) = sqrt(mean(diff(hrv_window_beat).^2));
    
    
    %% Frequency domain
    %Make the subwindow for Welch's Spectogram estimation, a tenth of the
    %length of the overall window
    sub_win = floor((window_end(wind_idx)-window_start(wind_idx))/10);
    if sub_win >1
        %Each section is windowed with a Hamming window
        [pxx, freq] = pwelch(detrend(hrv_window, 0),sub_win,[],[], fs_new);
        
        HRV_stats.VLF_power(wind_idx) = sum(pxx(and(freq> VLFP_low, freq< VLFP_high)));
        
        %get normalisation factor
        HRV_stats.T_power(wind_idx) = sum(pxx( freq< HFP_high));
        
        %Integrate the power in the frequency bands
        %     HRV_stats.nLF_power(wind_idx) = sum(pxx(and(freq> LFP_low, freq< LFP_high)))/HRV_stats.T_power(wind_idx);
        %     HRV_stats.nHF_power(wind_idx) = sum(pxx(and(freq> HFP_low, freq< HFP_high)))/HRV_stats.T_power(wind_idx);
        HRV_stats.nLF_power(wind_idx) = sum(pxx(and(freq> LFP_low, freq< LFP_high)));
        HRV_stats.nHF_power(wind_idx) = sum(pxx(and(freq> HFP_low, freq< HFP_high)));
    end
    %% Non linear
    
    %TODO Entropy
    
    %Do Poincare analysis
    X = [hrv_window_beat(1:end-1)';hrv_window_beat(2:end)'];
    if isempty(X)
        continue
    end
    X_rot = R * X;
    HRV_stats.poincare_SDD1(wind_idx) = std(X_rot(1, :));
    HRV_stats.poincare_SDD2(wind_idx) = std(X_rot(2, :));
    
    
    %     ellipse_t= fit_ellipse(hrv_window_beat(1:end-1), hrv_window_beat(2:end), fig);
end

HRV_stats.LF_div_HF  = HRV_stats.nHF_power./HRV_stats.nLF_power;


%should vary with these - but not exactly equal
HRV_stats.poincare_area = HRV_stats.poincare_SDD1 .* HRV_stats.poincare_SDD2;
HRV_stats.poincare_SDD1_SDD2 = HRV_stats.poincare_SDD1./ HRV_stats.poincare_SDD2; %Is correlated to LF/HF


%% PLOT
if configs.plot_flag
    
    fields = fieldnames(HRV_stats);
    fields_labels = strrep(fields, '_', ' ');
    window_centres = (t_window_start + t_window_end)/2;
    window_centres = window_centres/60;
    
    
    colours = func.aux_functions.define_colours();
    figure('Position', [3061         -66        1809         897])
    num_cols = 4;
    num_rows = ceil((1+length(fields))/num_cols);
    
    plot_idx = 1;
    
    ax(plot_idx) = subplot(num_rows, num_cols, plot_idx);
    plot_idx = plot_idx +1;
    plot(RR_times/60, RR_ts, '-ok')
    hold on
    plot(RR.t/60, RR.ts,'Color', colours.red)
    legend('Original', 'Interpolated')
    ylabel(configs.signal_name)
    xlabel('Time (mins)')
    
    hold on
    for p_idx = 1:length(fields)
        ax(p_idx+1) = subplot(num_rows, num_cols, p_idx +1);
        eval(['val = HRV_stats.',fields{p_idx},';'])
        plot(window_centres, val, '-o', 'Color', colours.green)
        ylabel(fields_labels{p_idx})
        xlabel('Time (mins)')
    end
    
    linkaxes(ax,'x')
    func.plot.tightfig();
end
end

