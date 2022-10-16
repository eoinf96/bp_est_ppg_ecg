% Function to get HRV features from a time series of RR intervals
% --
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
function [HRV_stats] = get_HRV_features(RR_ts, RR_times, window_start_times, window_end_times, opts)
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
    %     HRV_stats.samp_entropy = nan(1,1); % TODO!!
    HRV_stats.poincare_SDD1 = nan(1,1); %thought to be equivalent to RMSSD - short term HRV (sympathetic?)%Vasoconstriction is influenced by the low frequency sympathetic nervous system
    HRV_stats.poincare_SDD2 = nan(1,1);
    HRV_stats.poincare_area = nan(1,1);
    HRV_stats.poincare_SDD1_SDD2 = nan(1,1);
    return
end

if nargin < 3 || isempty(window_start_times) || isempty(window_end_times)
    %Recommended by the HR variability task force Camm et al 1996
    window_length_time = 60;
    window_step_time = 60;
    window_start_times = 0:window_step_time:(RR_times(end)-window_step_time);
    window_end_times = window_start_times + window_length_time;
end
% Set default params
if nargin < 5 || isempty(opts)
   opts = struct(); 
end
default_opts.plot_flag = 0;
default_opts.signal_name = 'RR interval';
opts = func.aux_functions.update_with_default_opts(opts, default_opts);


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
if isempty(window_start_times)
    window_start_times = RR_times(1);
    window_end_times = RR_times(end);
end

window_start = floor(window_start_times *fs_new)+1;
window_end = floor(window_end_times *fs_new);
window_end(window_end > length(RR.t)) = length(RR.t);
num_windows = length(window_start_times);


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
HRV_stats.samp_entropy = nan(num_windows,1); % TODO!!
HRV_stats.poincare_SDD1 = nan(num_windows,1); %thought to be equivalent to RMSSD - short term HRV (sympathetic?)
HRV_stats.poincare_SDD2 = nan(num_windows,1);
% poincare_area = nan(num_windows,1);

%Rotation matrix for use in poincare analysis
R = [cos(-pi/4) -sin(-pi/4);sin(-pi/4)  cos(-pi/4)] ;


for wind_idx = 1:num_windows
    hrv_window_beat = RR.ts_beat(and(RR.t_beat> window_start_times(wind_idx),RR.t_beat < window_end_times(wind_idx) ));
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
if opts.plot_flag
    
    fields = fieldnames(HRV_stats);
    fields_labels = strrep(fields, '_', ' ');
    window_centres = (window_start_times + window_end_times)/2;
    window_centres = window_centres/60;
    
    
    colours = constants_def('Colours');
    figure
    num_cols = 4;
    num_rows = ceil((1+length(fields))/num_cols);
    
    plot_idx = 1;
    
    ax(plot_idx) = subplot(num_rows, num_cols, plot_idx);
    plot_idx = plot_idx +1;
    plot(RR_times/60, RR_ts, '-ok')
    hold on
    plot(RR.t/60, RR.ts,'Color', colours.red)
    legend('Original', 'Interpolated')
    ylabel(opts.signal_name)
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
    set(gcf, 'Position', [3061         -66        1809         897])
    
end


% HRV_stats.t = (window_start_times + window_end_times)/2;;

end

