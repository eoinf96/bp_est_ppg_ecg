function PAT = pat_beat_ecg(PPG, ECG, plot_flag, fig_h, configs)
% pat_beat_ecg computed beat by beat pulse arrival time (PAT) from an
% input ECG and PPG waveform. PAT is computed as time from ECG to the first 
% PPG fiducial point before the next ECG. This function is typically called from
% get_PAT_beat which initialises the PPG.marker field
%
%
% INPUT: ECG: ECG time series
%        PPG: PPG time series
%        plot_flag - This saves a video of the computation of each PAT beat value. Useful for visualisation of
%        fig_h: Figure handle to plot detail on
%        configs: configs struct, see below for details
%
% OUTPUT: PAT: Beat by beat PAT time series. 
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
% ﻿Finnegan, E., Davidson, S., Harford, M., Jorge, J., Watkinson, P., Young, D., Tarassenko, L., & Villarroel, M. (2021). Pulse arrival time as a surrogate of blood pressure. Scientific Reports, 11(1), 1–21.
narginchk(2, inf)
if nargin < 3
    plot_flag = 0;
end
if nargin < 4
    fig_h = [];
end
if nargin < 5
    configs = struct();
end
if plot_flag && isempty(fig_h)
    fig_h = figure('visible', 'off');
end
default_configs.sqi_threshold =0.8; 
default_configs.file_save_loc = '.'; % Where you wanna save the video produced by plot_detail_flag
configs = func.aux_functions.update_with_default_opts(configs, default_configs);
%% Error checks
%Check ECG and PPG have been processed
if or(~isfield(PPG, 'peaks'), ~isfield(ECG, 'peaks'))
    error('Either the PPG or the ECG is not processed properly')
end

% Check for marker and assign default value
if ~isfield(PPG, 'marker')
   warning('No fiducial point included in the PPG struct. Will use onsets as default') 
   PPG.marker.ind = PPG.onsets;
   PPG.marker.t = PPG.t(PPG.marker.ind);
end
%% Set up video plotting
if plot_flag
    func.progress_bar.textprogressbar('Writing video '); 
    vid = VideoWriter([configs.file_save_loc, '/', PPG.record_name] );
    open(vid)
end
%% Initialise PAT struct
PAT.name = 'PAT';
PAT.units = 's';
if isfield(PPG, 'source'); PAT.source = ECG.source; end

num_of_beats = length(ECG.peaks) - 1;

% One PAT value for each R-wave
PAT.ts_beat = nan(num_of_beats, 1);
PAT.t_beat = ECG.t(ECG.peaks(1:num_of_beats));
PAT.sqi_beat = ones(num_of_beats, 1);
%% Go through each beat and extract PAT
for ecg_beat_index = 1:num_of_beats
    
    if ECG.sqi_beat(ecg_beat_index) < configs.sqi_threshold
        continue
    end
    
    ECG_peak_start = ECG.t(ECG.peaks(ecg_beat_index));
    ECG_peak_end= ECG.t(ECG.peaks(ecg_beat_index+1));
    
    %Is there a  PPG beat in between?
    ppg_beat_index = find(and(PPG.marker.t> ECG_peak_start,PPG.marker.t< ECG_peak_end));
    
    if isempty(ppg_beat_index) || length(ppg_beat_index) > 1
        PAT.sqi_beat(ecg_beat_index) = 0;
        continue
    end
    
    ppg_beat_index = ppg_beat_index(1);

    PAT.sqi_beat(ecg_beat_index) = min(ECG.sqi_beat(ecg_beat_index), PPG.sqi_beat(ppg_beat_index));
    PAT.ts_beat(ecg_beat_index) = PPG.marker.t(ppg_beat_index) - ECG_peak_start;
    %% Plot detial
    if plot_flag
        func.progress_bar.textprogressbar(100*ecg_beat_index/num_of_beats);
        
        %plot window and physiological range
        plot_beat(ECG, PPG, ecg_beat_index,ppg_beat_index, fig_h, PAT)
%         
        drawnow;
        frame = getframe(gcf);
        size(frame.cdata);
        writeVideo(vid, frame);
    end
    
end
if plot_flag
    func.progress_bar.textprogressbar('Done');
    close(vid)
end

%Remove nans where there are no beats within the range
loc_nan = isnan(PAT.ts_beat);
PAT.t_beat(loc_nan) = [];
PAT.ts_beat(loc_nan) = [];
PAT.sqi_beat(loc_nan) = [];

end



%% plot_detail
function plot_beat(ECG, PPG, ecg_beat_index,ppg_beat_index, fig_h, PAT)
%Inbuilt function to plot each window of PAT values
colours = constants_def('COLOUR');


qrs_duration = 0.25;%s
qrs_duration_samples = ceil(qrs_duration * ECG.fs);
half_qrs_duration_samples = ceil(qrs_duration_samples/2);

% Need to centre the ECG peak on zero
%%%%%%%%%%%%%%%

%first ECG beat
ECG_start = max(1,ECG.peaks(ecg_beat_index) - half_qrs_duration_samples);
%so we can centre the first ECG on 0
ECG_start_time = ECG.t(ECG.peaks(ecg_beat_index));
ECG_end = min(ECG.peaks(ecg_beat_index) + half_qrs_duration_samples, length(ECG.ts));

ECG_beat_f.ts = detrend(ECG.ts(ECG_start:ECG_end), 0);
%Align peak to be at 1
scaling_factor = max(ECG_beat_f.ts);
ECG_beat_f.ts = ECG_beat_f.ts/scaling_factor;


ECG_beat_f.t = ECG.t(ECG_start:ECG_end) - ECG_start_time;
ECG_beat_f.peaks = ECG.peaks(ecg_beat_index) - ECG_start+1;
%second ECG beat
ECG_start = max(0,ECG.peaks(ecg_beat_index +1) - half_qrs_duration_samples);
ECG_end = min(ECG.peaks(ecg_beat_index +1) + half_qrs_duration_samples, length(ECG.ts));

ECG_beat_s.ts = detrend(ECG.ts(ECG_start:ECG_end), 0);
ECG_beat_s.ts = ECG_beat_s.ts/scaling_factor;
ECG_beat_s.t = ECG.t(ECG_start:ECG_end) - ECG_start_time;

%%%%%%%%%%%%%%%

%Find the pleth we want to plot
beat_duration = 0.4;
max_time =  ECG.t(ECG_end);


t_start  = ECG_start_time ;
t_end = min(max_time + beat_duration, PPG.t(end));
index = find(and(PPG.t > t_start, PPG.t < t_end));

% PPG.ts = PPG.ts - PPG.ts(PPG.onsets(ppg_beat_index));
% PPG.ts = PPG.ts/PPG.ts(PPG.peaks(ppg_beat_index));

PPG_plot.ts = PPG.ts(index);
PPG_plot.ts = (PPG_plot.ts - min( PPG_plot.ts));
PPG_plot.ts = (PPG_plot.ts/ max( PPG_plot.ts));
PPG_plot.t = PPG.t(index) - ECG_start_time;

PPG_plot.markers = PPG.marker.ind(and(PPG.marker.ind > index(1), PPG.marker.ind <  index(end)));
PPG_plot.markers = PPG_plot.markers - index(1) +1;

PPG_plot.PAT_marker = PPG.marker.ind(ppg_beat_index);
PPG_plot.PAT_marker = PPG_plot.PAT_marker - index(1) +2;
%%

clf( fig_h );
%     figure( fig_h );

set(fig_h,'defaultAxesColorOrder',[colours.blue; colours.green]);


num_rows = 7;
num_cols = 4;

subplot(num_rows, num_cols, (1:num_cols))
%plot PAT
plot(PAT.t_beat, PAT.ts_beat, '-', 'Color', colours.blue, 'LineWidth', 1)
ylabel('PAT (s)')
xlim([PAT.t_beat(1), PAT.t_beat(end)])

subplot(num_rows, num_cols, (num_cols+1:2*num_cols))
plot(PAT.t_beat, PAT.ts_beat, '-*', 'Color', colours.blue, 'LineWidth', 1.3)
ylabel('PAT (s)')
half_window_time = 15;
t_win_start = max(PAT.t_beat(1), ECG_start_time - half_window_time);
t_win_end = min(PAT.t_beat(end), ECG_start_time + half_window_time);
xlim([t_win_start, t_win_end])


%set the xlim as the actual entire time span of the signal



subplot(num_rows, num_cols, (3*num_cols+1):(num_rows * num_cols))

%%
yyaxis left
plot(ECG_beat_f.t, ECG_beat_f.ts, 'Color', colours.blue,'LineWidth',3)
%     ylabel('ECG (mV)')
hold on
plot(ECG_beat_s.t, ECG_beat_s.ts, 'Color', colours.blue,'LineWidth',3, 'LineStyle', '-')
scatter(ECG_beat_f.t(ECG_beat_f.peaks), ECG_beat_f.ts(ECG_beat_f.peaks), 40, 'MarkerEdgeColor',colours.red,'MarkerFaceColor',colours.red)
ylim([min(min(ECG_beat_s.ts), min(ECG_beat_f.ts)), 1])
yticks([])

yyaxis right
plot(PPG_plot.t, PPG_plot.ts, 'Color', colours.green,'LineWidth',3)
scatter(PPG_plot.t(PPG_plot.markers), PPG_plot.ts(PPG_plot.markers), 40, 'MarkerEdgeColor',colours.purple,'MarkerFaceColor',colours.purple)

scatter(PPG_plot.t(PPG_plot.markers(1)), PPG_plot.ts(PPG_plot.PAT_marker), 40, 'MarkerEdgeColor',colours.yellow,'MarkerFaceColor',colours.yellow)

ylim([-0.1 1.1])
yticks([])


%     ylabel('PPG (adu)')

xlim([-qrs_duration/2, max(PPG_plot.t)])

set(gca, 'FontSize', 15)

yyaxis left
starting_x = 0;
ending_x = PPG_plot.t(PPG_plot.PAT_marker);
PAT = ending_x;
y = ECG_beat_f.ts(ECG_beat_f.peaks);
normalised1 = coords_to_pos(starting_x, y);
normalised2 = coords_to_pos(ending_x, y);
annotation('doublearrow',[normalised1(1) ,normalised2(1)],[normalised1(2), normalised2(2)]);


%     title(['t= ', num2str(ECG_start_time), ' (s)       SQI:', num2str(SQI)])
title(['PAT = ' num2str(PAT), ' (s)                  t= ', num2str(ECG_start_time), ' (s)'])
xticks([])
%     xlabel('Time (s)')

%     leg1.FontSize = 10;

end


function normalised  = coords_to_pos(x_coord, y_coord)

pos = get(gca, 'Position');
normalised(1) = (x_coord - min(xlim))/diff(xlim) * pos(3) + pos(1);
normalised(2) = (y_coord - min(ylim))/diff(ylim) * pos(4) + pos(2);

end