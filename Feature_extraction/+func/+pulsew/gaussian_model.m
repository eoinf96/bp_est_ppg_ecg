function [pts, rmse_error] = gaussian_model(ts,t,onsets,sqi_beat,config,do_plot)
% This function loops through each PPG beat and fits the sum of four
% gaussian model to each beat
% 
% INPUT:
%       ts -- PPG time series vecor
%       t -- time vector
%       onsets -- indices of onsets of each PPG pulse
%       sqi_beat -- sqi of each beat
%       config -- configs struct, see below for details
%       do_plot -- flag for plotting 
% OUTPUT:
%       pts -- struct containing amp, sigma and mu of each gaussian for each pulse
%       rmse_error -- Fitting error for each beat
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
% - Couceiro, R., Carvalho, P., Paiva, R.P., Henriques, J., Antunes, M., Quintal, I. and Mühlsteff, J., 2012, August. Multi-Gaussian fitting for the assessment of left ventricular ejection time from the Photoplethysmogram. In 2012 Annual International Conference of the IEEE Engineering in Medicine and Biology Society (pp. 3951-3954). IEEE.
%

if nargin < 5 || isempty(config)
    config = struct();
end
default_config.do_normalise = false; 
default_config.continue_points = true; %Whether we attempt to force consistency by using previous starting points as current starting points
default_config.error_threshold = 0.03; % Maximum allowable fitting error 
config = func.aux_functions.update_with_default_opts(config, default_config);

if nargin < 6
   do_plot = false; 
end
%% Set up
num_beats = length(onsets)-1;

names = {'g1', 'g2', 'g3', 'g4'};

for idx = 1:length(names)
    pts.(names{idx}).amp = nan(num_beats,1);
    pts.(names{idx}).mu = nan(num_beats,1);
    pts.(names{idx}).sigma = nan(num_beats,1);
end

% Define single Gaussian function
global gaussian
gaussian = @(b,x) b(1) * exp(-(x - b(2)).^2/b(3));


%initial parameters to start optimisation from
if config.do_normalise
    sigma_guess = 0.01;
else
    sigma_guess = median(diff(t(onsets)))/100;
end

standard_guess = [1/1., 1/5, sigma_guess, ...
        1/3, 2*1/5, sigma_guess,...
        1/4, 3*1/5, sigma_guess,...
        1/6, 4*1/5, sigma_guess];

       
%Store containing RMSE values of each beat
rmse_error = nan(num_beats,1);
%% Run loop
for pulse_no = 1 : num_beats
%     tic
    if sqi_beat(pulse_no) == 0
        continue
    end
    
    %% Get current pulse
    
    curr = [];
    %get current pulse
    curr_els = onsets(pulse_no):onsets(pulse_no+1);
    curr.ts = ts(curr_els);
    curr.t = t(curr_els);
    curr.t = curr.t - curr.t(1);
    
    % Correct for low frequency baseline drift in a single beat
    correction_line = linspace(curr.ts(1), curr.ts(end), length(curr.ts));
    curr.ts = curr.ts - correction_line' + curr.ts(1);
    
    %Normalise
    curr.ts_norm = curr.ts - min(curr.ts);
    curr.ts_norm = curr.ts_norm /max(curr.ts_norm);
    if config.do_normalise
        %Normalise time
        curr.t = curr.t - min(curr.t);
        curr.t = curr.t /max(curr.t);
    else
        t_range = curr.t(end);
        standard_guess = [1/1.1, 0.8*t_range/5, sigma_guess, ...
            1/3, 2*t_range/5, sigma_guess,...
            1/4, 3*t_range/5, sigma_guess,...
            1/6, 4*t_range/5, sigma_guess];
    end
    
    %% fit using previous opt as starting point
    if ~exist('p_opt', 'var')
        p_opt = standard_guess;
    end
    if config.continue_points
        use_standard = 0;
        try
            p_opt_prev = local_fit_model(p_opt,curr.t, curr.ts_norm);
        catch ME
            %Error in fit
            fprintf(sprintf('Error: %s \n',ME.message ))
            use_standard =1;
            p_opt_prev = p_opt;
        end
    else
        use_standard = 1;
    end
    %% Fit using standard start
    try
        p_opt_stand = local_fit_model(standard_guess,curr.t, curr.ts_norm);
    catch ME
        %Error in fit
        fprintf(sprintf('Error: %s \n',ME.message ))
        if use_standard
            continue
        end
        p_opt_stand = standard_guess;
    end
    %% Compare the guesses using previous and standard guess
    if use_standard ==1
        p_opt = p_opt_stand;
        rmse_error(pulse_no)=  get_error(p_opt, curr.t, curr.ts_norm); 
    else
        %Compare RMSE values
        [RMSE_prev] = get_error(p_opt_prev, curr.t, curr.ts_norm);
        [RMSE_stand] = get_error(p_opt_stand, curr.t, curr.ts_norm);
        
        if RMSE_prev < RMSE_stand
            p_opt = p_opt_prev;
            rmse_error(pulse_no)=  RMSE_prev; 
        else
            p_opt = p_opt_stand;
            rmse_error(pulse_no)=  RMSE_stand; 
        end
    end
    
    
    %% Assign new points
    for idx = 1:length(names)
        pts.(names{idx}).amp(pulse_no)      = p_opt(3*(idx-1)+1);
        pts.(names{idx}).mu(pulse_no)       = p_opt(3*(idx-1)+2);
        pts.(names{idx}).sigma(pulse_no)    = p_opt(3*(idx-1)+3);
    end

    
    for idx = 1:length(names)
        eval(['pts.',names{idx},'.amp(pulse_no) = p_opt(',num2str(3*(idx-1)+1),');'])
        eval(['pts.',names{idx},'.mu(pulse_no) = p_opt(',num2str(3*(idx-1)+2),');'])
        eval(['pts.',names{idx},'.sigma(pulse_no) = p_opt(',num2str(3*(idx-1)+3),');'])
    end
    
    if do_plot
        figure
        g1 = gaussian(p_opt(1:3), curr.t);
        g2 = gaussian(p_opt(4:6), curr.t);
        g3 = gaussian(p_opt(7:9), curr.t);
        g4 = gaussian(p_opt(10:12), curr.t);
        overall = g1 + g2 + g3 + g4;
        
        plot(curr.t, curr.ts_norm)
        hold on
        plot(curr.t, g1, 'r-', 'LineWidth', 2);
        plot(curr.t, g2, 'b-', 'LineWidth', 2);
        plot(curr.t, g3, 'y-', 'LineWidth', 2);
        plot(curr.t, g4, 'm-', 'LineWidth', 2);
        plot(curr.t, overall, 'k--', 'LineWidth', 3);
    end
    close
%     toc
end
%% Error handling
% Any Gaussian fit with an error greater than the threshold is set to nan;
set_rmse_gauss_func = @(x) set_rmse_gauss(x, rmse_error, config.error_threshold);
pts = structfun(set_rmse_gauss_func, pts, 'UniformOutput', false);
end

%% Local functions
function p_opt = local_fit_model(p0,t, ts)
% This function runs the Gaussian model fitting using initial parameters p0
func_in_line = @(p, x) get_overall(p, x);
lb = zeros(1,12);
ub = repmat([1, 2, inf], 1,4);
options = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt', 'Display','off');%, 'UseParallel', true);
p_opt = lsqcurvefit(func_in_line,p0, t,ts,lb,ub,options);
end



function out = get_overall( b,t)
%This function sets up the objective function for the optimoptions to optimise
global gaussian 
%Check that they are in the correct order
mu_1 = b(2);
mu_2 = b(5);
mu_3 = b(8);
mu_4 = b(11);
if mu_1 < mu_2 && mu_2 < mu_3 && mu_3 < mu_4
    modelfun = @(b,x) gaussian(b(1:3), x(:,1)) + gaussian(b(4:6), x(:,1))+ gaussian(b(7:9), x(:,1))+ gaussian(b(10:12), x(:,1));
    out = modelfun(b, t);
else
    out = 1*10^(28)* ones(size(t)); % Big number
end
end


function [RMSE]= get_error(b, t, ts)
%This function computes the error of the fit
est_fit = get_overall(b,t);
RMSE = sqrt(mean((est_fit - ts).^2));
% NRMSE = RMSE./range(ts);
end

function out = set_rmse_gauss(in, rmse, error_threshold)
   out= in;
   out.amp(rmse > error_threshold) = nan;
   out.mu(rmse > error_threshold) = nan;
   out.sigma(rmse > error_threshold) = nan;
end

