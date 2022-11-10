function derivs = get_derivs(ts, fs, configs)
% This function locates the characteristic fiducial points of the PPG
% Code adapted from P.Charlton:
% https://github.com/peterhcharlton/pulse-analyse
%
% INPUT:ts: signal
%       fs: sampling frequency 
%       configs: configs struct, see below for details
% 
% OUTPUT: pts: struct of derivatives of signal
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
% - Abraham Savitzky and Marcel J E Golay. “Smoothing and differentiation of data by simplified least squares procedures.” In: Analytical chemistry 36.8 (1964), pp. 1627–1639.

narginchk(2, inf);
if nargin < 3
    configs = struct();
end
default_configs.s_g_filter_len = 7;
default_configs.max_deriv_order = 4; % Maximum derivative order to compute -- int up to 8
configs = func.aux_functions.update_with_default_opts(configs, default_configs);
deriv_names = {'first', 'second', 'third', 'fourth', 'fifth', 'sixth', 'seventh', 'eighth'};
%% Error check
if configs.max_deriv_order > length(deriv_names)
    warning('Max deriv order too high set to the default value: %.0f \n', length(deriv_names))
    configs.max_deriv_order = length(deriv_names);
end
%% Get derivs
dt = 1/fs;
for der_i = 1:configs.max_deriv_order
    derivs.(deriv_names{der_i})  = func.waveform.savitzky_golay_deriv(ts, der_i, configs.s_g_filter_len)./(dt.^(der_i));
end
end