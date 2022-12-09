function PAT = get_pat_beat(ECG,PPG, distal_point, configs, plot_flag)
%  This function computes beat by beat PAT estimates as the time delay
%  between fiducial points on the ECG and PPG signals. 
%
% INPUT: ECG: ECG time series
%        PPG: PPG time series
%        distal_point: name of PPG distal point to use - default is 'tangent'
%        configs: configs struct, see below for details
%        plot_flag
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
%% Proces inputs
narginchk(2, inf);
if nargin < 3 || isempty(distal_point)
    distal_point = 'tangent';% tangent has been shown to be the least affected by wave reflections
end
if nargin < 4 || isempty(configs)
    configs = struct();
end
if nargin < 5 || isempty(plot_flag)
    plot_flag = false;
end
% Set the default configs
default_configs.PAT_range_flag =0; %This indicates whether to compute PAT between 2 ecg pulses (0) or witihin a defined range
default_configs.PAT_min_value =0.3; % Only required if PAT_range_flag = 1
default_configs.PAT_max_value =0.8; % Only required if PAT_range_flag = 1
default_configs.sqi_threshold =0.8; 
default_configs.plot_detail_flag = 0; % Little bit of extra plotting fun
default_configs.file_save_loc = '.'; % Where you wanna save the video produced by plot_detail_flag

% Update configs
configs = func.aux_functions.update_with_default_opts(configs, default_configs);
%% Error checks
%Check ECG and PPG have been processed
if or(~isfield(PPG, 'peaks'), ~isfield(ECG, 'peaks'))
    error('Either the PPG or the ECG is not processed properly')
end
%% Set variables for PPG fiducial point
%We follow Elgendi's descirption of fiducial points:
%Mohamed Elgendi, Yongbo Liang, and Rabab Ward. Toward Generating More Diagnostic Features from Photoplethysmogram Waveforms. Diseases, 6(1):20, 2018. ISSN 2079-9721. doi: 10.3390/ diseases6010020.
switch distal_point
    case 'peaks'
%         PPG.marker = PPG.fid_pts.s;
        PPG.marker.ind = PPG.peaks;
        PPG.marker.t = PPG.t(PPG.peaks);
        PPG.marker.amp = PPG.ts(PPG.peaks);
    case 'onsets'
        PPG.marker.amp = PPG.fid_pts.f1.amp;
        PPG.marker.amp(end+1) = PPG.fid_pts.f2.amp(end);
        PPG.marker.ind = PPG.fid_pts.f1.ind;
        PPG.marker.ind(end+1) = PPG.fid_pts.f2.ind(end);
        PPG.marker.t = PPG.fid_pts.f1.t;
        PPG.marker.t(end+1) = PPG.fid_pts.f2.t(end);
        PPG.sqi_beat(end+1) = PPG.sqi_beat(end);
    case 'halfpoint'
        PPG.marker = PPG.fid_pts.halfpoint;
    case 'tangent'
        PPG.marker = PPG.fid_pts.tangent;
        PPG.marker.ind = round(PPG.marker.ind);
    case 'max_deriv'
        PPG.marker = PPG.fid_pts.W;
        PPG.marker.amp = PPG.ts(PPG.marker.ind(~isnan(PPG.marker.ind)));
    case 'max_sec_deriv'
        PPG.marker = PPG.fid_pts.a;
        PPG.marker.amp = PPG.ts(PPG.marker.ind(~isnan(PPG.marker.ind)));
    otherwise
        error('Unkown distal point')
end

PPG.marker.t(PPG.sqi_beat ==0) = nan;
remove_marker = isnan(PPG.marker.t);

PPG.sqi_beat(remove_marker) = [];
PPG.marker.ind(remove_marker) = [];
PPG.marker.t(remove_marker) = [];
PPG.marker.amp(remove_marker) = [];
%% Get PAT



if configs.PAT_range_flag
    PAT = func.pat.pat_beat_range(PPG, ECG,configs.PAT_min_value, configs.PAT_max_value ,configs.plot_detail_flag, [], configs );
else
    PAT = func.pat.pat_beat_ecg(PPG, ECG,configs.plot_detail_flag, [] );
end

PAT.distal_point = distal_point;
if isfield(ECG, 'start_datetime') || isfield(PPG, 'start_datetime')
    PAT.start_datetime = ECG.start_datetime;
    PAT.dt_beat = PAT.start_datetime + seconds(PAT.t_beat);
end
% PAT= func.aux_functions.update_good(PAT);
%% plot beat by beat PAT values
if plot_flag
   
    do_mins = ECG.t(end) > 10*60;
    if do_mins % Plot mins or seconds
        factor = 60;
    else
        factor = 1;
    end
    
    colours = func.aux_functions.define_colours;
    %%
    figure('Position', [377   163   862   776])
    num_rows = 11;
    index = 1;
    ax_index = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ax(ax_index) = subplot(num_rows,1,index:index+1);
    ax_index= ax_index +1;
    index = index+2;
    plot(ECG.t/factor, ECG.ts, 'LineWidth',1.2, 'Color', colours.black)
    hold on
    h1 = scatter(ECG.t(ECG.peaks)/factor, ECG.ts(ECG.peaks), '*','MarkerEdgeColor', colours.red);
    ylabel('ECG (mV)');
    legend(h1, 'R-peaks')
    set(gca,'xticklabel',{[]})
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax(ax_index) = subplot(num_rows,1,index);
    ax_index= ax_index +1;
    index = index+1;
    plot(ECG.t/factor, ECG.sqi, 'LineWidth',1.2,'color',colours.black)
    hold on
    line([ECG.t(1)/factor, ECG.t(end)/factor], [configs.sqi_threshold, configs.sqi_threshold],'LineStyle','--', 'color', colours.yellow)
    ylabel('ECG sqi');
    
    ylim([0 1.1])
    set(gca,'xticklabel',{[]})

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    ax(ax_index) = subplot(num_rows,1,index:index+1);
    ax_index= ax_index +1;
    index = index+2;
    plot(PPG.t/factor, PPG.ts, 'LineWidth',1.2, 'Color', colours.black)
    hold on
    h1 = scatter(PPG.marker.t/factor, PPG.marker.amp, '*','MarkerEdgeColor', colours.red);
    ylabel('PPG (a.d.u)');
    legend(h1,  ['Fiducial point - ', distal_point])
    set(gca,'xticklabel',{[]})

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    ax(ax_index) = subplot(num_rows,1,index);
    ax_index= ax_index +1;
    index = index+1;
    plot(PPG.t/factor, PPG.sqi, 'LineWidth',1.2,'Color',colours.black)
    hold on
    line([PPG.t(1)/factor, PPG.t(end)/factor], [configs.sqi_threshold, configs.sqi_threshold],'LineStyle','--', 'color', colours.yellow)
    ylabel('PPG sqi');
    ylim([0 1.1])
    set(gca,'xticklabel',{[]})


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    ax(ax_index) = subplot(num_rows,1,index:index+1);
    ax_index= ax_index +1;
    index = index+2;
    plot(PAT.t_beat/factor, PAT.ts_beat, '-*', 'LineWidth',0.8, 'Color', colours.black)
    if configs.PAT_range_flag
        hold on
        line([PAT.t_beat(1)/factor,PAT.t_beat(end)/factor], [configs.PAT_min_value, configs.PAT_min_value], 'Color', colours.yellow, 'LineStyle', '--')
        line([PAT.t_beat(1)/factor,PAT.t_beat(end)/factor], [configs.PAT_max_value, configs.PAT_max_value], 'Color', colours.yellow, 'LineStyle', '--')        
        ylim([configs.PAT_min_value,configs.PAT_max_value])
        legend('Beat by beat PAT', 'Max PAT', 'Min PAT');
    end
    ylabel('PAT');
    set(gca,'xticklabel',{[]})
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    ax(ax_index) = subplot(num_rows,1,index);
    ax_index= ax_index +1;
    index = index+1;
    plot(PAT.t_beat/factor, PAT.sqi_beat, 'LineWidth',1.2, 'LineWidth',1.2,'Color',colours.black)
    hold on
    line([PAT.t_beat(1)/factor, PAT.t_beat(end)/factor], [configs.sqi_threshold, configs.sqi_threshold],'LineStyle','--', 'color', colours.yellow)
    ylabel('PAT sqi');
    ylim([0 1.1])
    set(gca,'xticklabel',{[]})


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ax(ax_index) = subplot(num_rows,1,index:index+1);
    if ~isfield(PAT, 'good_beat'); PAT.good_beat = PAT.ts_beat; PAT.good_beat(PAT.sqi_beat == 0 ) = nan;     end
    plot(PAT.t_beat/factor, PAT.good_beat, '-*', 'LineWidth',0.8, 'Color', colours.black)
    ylabel('PAT good');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    linkaxes(ax, 'x')
    
    if do_mins
        xlabel('Time (mins)');
    else
        xlabel('Time (secs)');
    end

    set(findall(gcf,'-property','FontSize'),'FontSize',12)
    func.plot.tightfig;
end

end

