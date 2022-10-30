function configs = return_configs()
% This function returns configs struct
% 
% OUTPUT: pts: struct defining the locations, timings and amplitudes of all fiducial points
%         norm_pts: struct defining the locations, timings and amplitudes of all fiducial points using normalised PPG pulses -- only returned if config.do_normalise == 1
%           derivs: struct of PPG derivatives
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


    %%% Fid point config
    configs.PPG.fid_point.do_fid = 1;
    configs.PPG.fid_point.gauss_continue_points = true;
    configs.PPG.fid_point.do_e_for_dic = true;
    configs.PPG.fid_point.do_normalise = 1;

    %change as necessary
%     configs.PPG.fid_point.fid_pt_names = {'a', 'b', 'c', 'd', 'e', 'f', 's', 'dia', 'dic', ...
%             'W', 'f1', 'f2', 'halfpoint', 'tangent', 'gauss'};  
    configs.PPG.fid_point.fid_pt_names = {'a', 'b', 'c', 'd', 'e', 'f', 's', 'dia', 'dic', ...
        'W', 'f1', 'f2', 'halfpoint', 'tangent'};  
        
end

