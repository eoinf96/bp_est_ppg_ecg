function PAT = pat_beat_range(PPG,...
                                ECG,...
                                PAT_min_value,...
                                PAT_max_value,...
                                plot_flag,...
                                fig_h,...
                                configs ...
                                )
% PAT_beat_range computed beat by beat pulse arrival time (PAT) from an
% input ECG and PPG waveform. PAT is computed as time from ECG to some PPG 
% point within a defined range. This function is typically called from
% get_PAT_beat which initialises the PPG.marker field

%% Proces inputs
narginchk(2, inf)
if nargin < 3
   PAT_min_value = 0;
end
if nargin < 4
    PAT_min_value = 0.5;
end
if nargin < 5
   plot_flag = 0; 
end
if nargin < 6
    fig_h = [];
end
if nargin < 7
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
    external.progress_bar.textprogressbar('Writing video '); 
    vid = VideoWriter([configs.file_save_loc, '/', PPG.record_name] );
    open(vid)
end
%% Initialise PAT struct
PAT.name = 'PAT';
PAT.units = 's';
PAT.source = ECG.source;

num_of_beats = length(ECG.peaks) - 1;

% One PAT value for each R-wave
PAT.ts_beat = nan(num_of_beats, 1);
PAT.t_beat = ECG.t(ECG.peaks(1:num_of_beats));
PAT.sqi_beat = ones(num_of_beats, 1);

marker_times = PPG.marker.t;
%% Go through each ECG beat and extract PAT
for ecg_beat_index = 1:num_of_beats
    
    if ECG.sqi_beat(ecg_beat_index) < configs.sqi_threshold
        continue
    end
    
    ECG_peak_start = PAT.t_beat(ecg_beat_index);
    
    %Choose first beat within the range
    ppg_beat_index = find(and(marker_times > ECG_peak_start +PAT_min_value, marker_times < ECG_peak_start +PAT_max_value), 1, 'first');
    %If there are no PAT values within this range then keep PAT as nan 
    %-- we will remove nan beats later
    if isempty(ppg_beat_index)
        PAT.sqi_beat(ecg_beat_index) = 0;
        continue
    end
    
    PAT.ts_beat(ecg_beat_index) = marker_times(ppg_beat_index) - ECG_peak_start;
    PAT.sqi_beat(ecg_beat_index) = min(ECG.sqi_beat(ecg_beat_index), PPG.sqi_beat(ppg_beat_index));
    

    if plot_flag
        external.progress_bar.textprogressbar(100*ecg_beat_index/num_of_beats);
        %plot window and physiological range
        plot_beat(ECG, PPG, ecg_beat_index,ppg_beat_index, PAT_max_value, PAT_min_value, fig_h, PAT.sqi_beat(ecg_beat_index))
        drawnow;
        frame = getframe(gcf);
        %             pause(0.05)
        size(frame.cdata);
        writeVideo(vid, frame);
    end
    
    %Remove the ppg beat from vec so that it cannot be used again
    marker_times(ppg_beat_index) = nan;
    
end
if plot_flag
    external.progress_bar.textprogressbar('Done');
    close(vid)
end



%Remove nans where there are no beats within the range
loc_nan = isnan(PAT.ts_beat);
PAT.t_beat(loc_nan) = [];
PAT.ts_beat(loc_nan) = [];
PAT.sqi_beat(loc_nan) = [];

end


%% plot_detail
function plot_beat(ECG, PPG, ecg_beat_index,ppg_beat_index, PAT_max_value, PAT_min_value, fig_h, SQI)
%Inbuilt function to plot each window of PAT values
colours = constants_def('COLOUR');

qrs_duration = 0.25;%s
qrs_duration_samples = ceil(qrs_duration * ECG.fs);
half_qrs_duration_samples = ceil(qrs_duration_samples/2);


ECG_start = max(0,ECG.peaks(ecg_beat_index) - half_qrs_duration_samples);
%so we can centre the first ECG on 0
ECG_start_time = ECG.t(ECG.peaks(ecg_beat_index));
ECG_end = min(ECG.peaks(ecg_beat_index) + half_qrs_duration_samples, length(ECG.ts));

ECG_beat.ts = detrend(ECG.ts(ECG_start:ECG_end), 0);
%Align peak to be at 1
scaling_factor = max(ECG_beat.ts);
ECG_beat.ts = ECG_beat.ts/scaling_factor;

ECG_beat.t = ECG.t(ECG_start:ECG_end) - ECG_start_time;
ECG_beat.peaks = ECG.peaks(ecg_beat_index) - ECG_start+1;


%Find the pleth we want to plot
beat_duration = 0.4;

t_start  = ECG_start_time ;
t_end = min(ECG_start_time + PAT_max_value + beat_duration, PPG.t(end));
index = find(and(PPG.t > t_start, PPG.t < t_end));


PPG_plot.ts = PPG.ts(index);
PPG_plot.t = PPG.t(index)- ECG_start_time;
% All the possible markers
PPG_plot.markers = PPG.marker.ind(and(PPG.marker.ind > index(1), PPG.marker.ind < index(end)));
PPG_plot.markers = PPG_plot.markers - index(1) +1;

% The marker that was chosen
PPG_plot.PAT_marker = PPG.marker.ind(ppg_beat_index);
PPG_plot.PAT_marker = PPG_plot.PAT_marker - index(1) +1;
%%

clf( fig_h );
%     figure( fig_h , 'visible', 'off');

set(fig_h,'defaultAxesColorOrder',[colours.blue; colours.green]);

%%
yyaxis left
plot(ECG_beat.t, ECG_beat.ts, 'Color', colours.blue,'LineWidth',3)
%     ylabel('ECG (mV)')
yticks([])
hold on
scatter(ECG_beat.t(ECG_beat.peaks), ECG_beat.ts(ECG_beat.peaks), 40, 'MarkerEdgeColor',colours.red,'MarkerFaceColor',colours.red)
yyaxis right
plot(PPG_plot.t, PPG_plot.ts, 'Color', colours.green,'LineWidth',3)
scatter(PPG_plot.t(PPG_plot.markers), PPG_plot.ts(PPG_plot.markers), 40, 'MarkerEdgeColor',colours.purple,'MarkerFaceColor',colours.purple)
scatter(PPG_plot.t(PPG_plot.PAT_marker), PPG_plot.ts(PPG_plot.PAT_marker), 40, 'MarkerEdgeColor',colours.yellow,'MarkerFaceColor',colours.yellow)

yticks([])
%     ylabel('PPG (adu)')

line([PAT_min_value PAT_min_value], [min(PPG_plot.ts), max(PPG_plot.ts)], 'Color', colours.yellow, 'LineStyle', '--','LineWidth',3)
line([PAT_max_value PAT_max_value], [min(PPG_plot.ts), max(PPG_plot.ts)], 'Color', colours.yellow, 'LineStyle', '--','LineWidth',3)
%     leg1 = legend('ECG', 'R peak', 'PPG','PPG markers' ,'Physiological range');
title(['t= ', num2str(ECG_start_time), ' (s)       SQI:', num2str(SQI)])

xlabel('Time (s)')
set(gca, 'FontSize', 15)


yyaxis left
starting_x = 0;
ending_x = PPG_plot.t(PPG_plot.PAT_marker);
y = ECG_beat.ts(ECG_beat.peaks);
normalised1 = coords_to_pos(starting_x, y);
normalised2 = coords_to_pos(ending_x, y);
annotation('doublearrow',[normalised1(1) ,normalised2(1)],[normalised1(2), normalised2(2)]);

end

function normalised  = coords_to_pos(x_coord, y_coord)

pos = get(gca, 'Position');
normalised(1) = (x_coord - min(xlim))/diff(xlim) * pos(3) + pos(1);
normalised(2) = (y_coord - min(ylim))/diff(ylim) * pos(4) + pos(2);

end