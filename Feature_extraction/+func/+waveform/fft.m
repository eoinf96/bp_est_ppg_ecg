function [pw, ff] = fft(ts, fs, configs, plot_flag)
% This function computes the Fast Fourier Transform of a signal 
%
% INPUT: ts: Signal in
%        fs: Sampling frequency of signal in
%        configs: configs struct, see below for details
%        plot_flag
%
% OUTPUT: pw: Power of the FFT
%         ff: Corresponding frequency
%
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
narginchk(2, inf)
%% Update opts
if nargin < 3 || isempty(configs)
    configs  = struct();
end
if nargin < 4
    plot_flag  = false;
end
default_configs.window_data = 0; % Whether to multiple the signal by a Hamming window in order to reduce the effect of discontinuities
default_configs.detrend =1;
configs = func.aux_functions.update_with_default_opts(configs, default_configs);
%% Preprocess
ts = ts(:);	
N = length(ts);
if configs.detrend 
    ts =detrend(ts, 0, 'omitnan');
end
if configs.window_data
    window_data = hamming( N );
    ts = ts .* window_data;
end
%% Compute FFT
N_FFT = 2^nextpow2( N );  
y = fft(ts, N_FFT);

% frequency vector
ff = (0:N_FFT/2)' * fs/N_FFT; 

% Take the magnitude of fft and scale it
pw_mag = abs( y(1:(N_FFT/2+1)) / N_FFT );
sq_mag = pw_mag .^2;
pw = sq_mag / norm(sq_mag);

%% Plot flag

if plot_flag
    figure('Position', [446   670   817   250]);
    subplot(1, 4, 1:3);
    plot(linspace(0, (length(ts)-1)/fs, length(ts)), ts, 'k')
    xlabel('t (secs)')
    ylabel('Sig')
    
    % Plot FFT
    subplot(1, 4, 4);
    plot(ff, pw, 'k')
    xlabel('f (Hz)')
    ylabel('FFT')
    
    func.plot.tightfig();
end
end

