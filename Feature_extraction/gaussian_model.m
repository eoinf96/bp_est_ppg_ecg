% Functions to fit the sum of four Gaussians to each PPG beat
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
%%
function [pts] = gaussian_model(ts,...
                                t,...
                                onsets,...
                                sqi_beat,...
                                params,...
                                do_plot)
% This function loops through each PPG beat and fits the sum of four
% gaussian model to each beat
% 
% Inputs:
%       ts -- PPG time series vecor
%       t -- time vector
%       onsets -- indices of onsets of each PPG pulse
%       sqi_beat -- sqi of each beat
%       do_normalise -- flag indicating whether to normalise the width of
%                       each beat to be 1
%       do_plot -- flag for plotting 
% Outputs:
%       pts -- struct containing amp, sigma and mu of each gaussian for
%       each pulse
if nargin < 5 || isempty(params)
    params = struct();
end
default_params.do_normalise = false; 
default_params.continue_points = true; %Whether we attempt to force consistency by using previous starting points
params = func.aux_functions.update_with_default_opts(params, default_params);

if nargin < 6
   do_plot = false; 
end
%% Set up
num_beats = length(onsets)-1;

names = {'g1', 'g2', 'g3', 'g4'};

for idx = 1:length(names)
    eval(['pts.',names{idx},'.amp = nan(num_beats,1);'])
    eval(['pts.',names{idx},'.mu = nan(num_beats,1);'])
    eval(['pts.',names{idx},'.sigma = nan(num_beats,1);'])
end

global gaussian
gaussian = @(b,x) b(1) * exp(-(x - b(2)).^2/b(3));


%initial parameters to start optimisation from
if params.do_normalise
    sigma_guess = 0.01;
else
    sigma_guess = median(diff(t(onsets)))/100;
end

standard_guess = [1/1., 1/5, sigma_guess, ...
        1/3, 2*1/5, sigma_guess,...
        1/4, 3*1/5, sigma_guess,...
        1/6, 4*1/5, sigma_guess];

       
%Store containing RMSE and normalised RMSE values of each beat
rmse_error = nan(num_beats,1);
%% Run loop
for pulse_no = 1 : num_beats
%     tic
    if sqi_beat(pulse_no) < 0.8
        continue
    end
    
    %% Get pulse
    
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
    if params.do_normalise
        %Normalise time
        curr.t = curr.t - min(curr.t);
        curr.t = curr.t /max(curr.t);
    end
    
    
    %Initial guess at params
    if ~params.do_normalise
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
    if params.continue_points
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

error_threshold = 0.03; % have slightly increased the error threshold to not be so strict
set_rmse_gauss_func = @(x) set_rmse_gauss(x, rmse_error, error_threshold);
pts = structfun(set_rmse_gauss_func, pts, 'UniformOutput', false);

pts.gauss_error = rmse_error;



end

%% Local functions
function p_opt = local_fit_model(p0,t, ts)
% This function runs the model fitting
func_in_line = @(p, x) get_overall(p, x);
lb = zeros(1,12);
ub = repmat([1, 2, inf], 1,4);
options = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt', 'Display','off');%, 'UseParallel', true);
p_opt = lsqcurvefit(func_in_line,p0, t,ts,lb,ub,options);
end



function out = get_overall( b,t)
%This function sets up the objective function for the optimoptions to
%optimise
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

