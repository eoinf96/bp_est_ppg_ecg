function [colours] = define_colours
% This function returns RGB colours typically used for plotting
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


% Define colours for plotting
    colours.blue        = [0.3010, 0.7450, 0.9330];
    colours.red         = [0.6350, 0.0780, 0.1840];
    colours.yellow      = [0.9290, 0.6940, 0.1250];
    colours.green       = [0.1, 0.7, 0.6];
    colours.purple      = [0.4940, 0.1840, 0.5560];
    colours.dark_green  = [0,0.3922,0];
%     colours.orange      = [1,0.6471,0];
    colours.orange      = [0.93, 0.53, 0.18];
    colours.dark_blue   = [0,0,0.8039];
    colours.bright_blue = [0, 0.75, 1];
    colours.black       = [0, 0, 0];
    colours.grey        = [0.4, 0.4, 0.4];
end

