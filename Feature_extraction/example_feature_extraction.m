%% Example code of feature extraction from the ECG and PPG
% Please cite the following when using this code base:
% Finnegan, E., Davidson, S., Harford, M., Jorge, J., Watkinson, Tarassenko, L., \& Villarroel, M. 
% Features from the photoplethysmogram and electrocardiogram for estimating changes in blood pressure.
% Submitted to Scientific reports
close all
clear all
clc
%% Load data
load('example_data/ECG.mat', 'ECG');
load('example_data/PPG_all.mat', 'PPG');

% The ECG and PPG have already been segmented and their sqi is already determined. The code for performing this is propietary at the moment. 
% All steps for segmenting the ECG and PPG are outlined in the associated
% publication. We note below some additionaly relevant publications:
% - Mauricio Villarroel et al. “Non-contact vital-sign monitoring of patients undergoing haemodialysis treatment”. In: Scientific reports 10.1 (2020), pp. 1–21.
% - Q. Li and G. D. Clifford. “Dynamic time warping and machine learning for signal quality assessment of pulsatile signals”. In: Physiological Measurement 33.9 (2012), pp. 1491–1501.
% - Q. Li, R. G. Mark, and G. D. Clifford. “Robust heart rate estimation from multiple asynchronous noisy sources using signal quality indices and a Kalman filter”. In: Physiological Measurement 29.1 (Jan. 2008), pp. 15–32.
%% Flags
flags.plot_ppg_ecg_overview = 0;
flags.plot_ppg_features_overview = 0;
%% Parms
params.window_size = 30;
params.window_step = 30;
%% Options
configs_ECG.do_HRV_feats = 1;
%% Visualise ECG and PPG

if flags.plot_ppg_ecg_overview 
colours = func.aux_functions.define_colours();
lw = 1.6;
figure('position', [84         508        1484         439])
num_rows = 4; num_cols = 1; p_idx = 1;
ax(p_idx) = subplot(num_rows, num_cols, p_idx);hold on; p_idx = p_idx+1;
plot(PPG.t/60, PPG.ts, 'Color', colours.black, 'LineWidth', lw);
xlabel('Time (mins)'); ylabel('PPG');
ax(p_idx) = subplot(num_rows, num_cols, p_idx);hold on; p_idx = p_idx+1;
plot(PPG.t(PPG.peaks)/60, PPG.sqi_beat, 'Color', colours.black, 'LineWidth', lw);
ylim([0, 1])
xlabel('Time (mins)'); ylabel('SQI_{PPG}');
ax(p_idx) = subplot(num_rows, num_cols, p_idx);hold on; p_idx = p_idx+1;
plot(ECG.t/60, ECG.ts, 'Color', colours.black, 'LineWidth', lw);
xlabel('Time (mins)'); ylabel('ECG');
ax(p_idx) = subplot(num_rows, num_cols, p_idx);hold on; p_idx = p_idx+1;
plot(ECG.t(ECG.peaks)/60, ECG.sqi_beat, 'Color', colours.black, 'LineWidth', lw);
xlabel('Time (mins)'); ylabel('SQI_{ECG}');
ylim([0, 1])
linkaxes(ax, 'x')
func.plot.tightfig;
end
%% Get features from PPG
%Load configs
configs = return_configs;
[PPG.fid_pts, PPG.norm_fid_pts, PPG.derivs]     = func.pulsew.get_ppg_fid_pts(PPG, configs.PPG.fid_point); 
PPG.pw_inds                                     = func.pulsew.get_ppg_indices(PPG, 0);
PPG.norm_pw_inds                                = func.pulsew.get_ppg_indices(PPG, 1);
%% Get PPG feature vector
t_window_start = 0:params.window_size:(PPG.t(end)-30);
t_window_end = t_window_start + params.window_step;
PPG_feature_table = func.pulsew.make_PPG_feature_vector(PPG, t_window_start, t_window_end);
%% Plot PPG featuress
if flags.plot_ppg_features_overview; func.pulsew.plots.plot_indices_and_gauss(PPG, 15); end
%% Get ECG feature vector
t_window_start = 0:params.window_size:(ECG.t(end)-30);
t_window_end = t_window_start + params.window_step;
ECG_feature_table = func.ecg.make_ECG_feature_vector(ECG, t_window_start, t_window_end, configs_ECG);
%% Get PAT
PAT = func.pat.get_pat_beat(ECG, PPG, 'tangent', [], 1);
PAT_feature_table = func.pat.make_pat_feature_vector(PAT, t_window_start, t_window_end);
%% Get total feature vector
feature_vector = [PPG_feature_table, ECG_feature_table, PAT_feature_table];
