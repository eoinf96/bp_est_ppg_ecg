function opts = update_with_default_opts(opts,default_opts)
% This function  updates any function configs with the default values
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
opt_names = fieldnames(default_opts); 
for name_idx = 1:length(opt_names)
    if ~isfield(opts,opt_names{name_idx} )
        opts.(opt_names{name_idx}) = default_opts.(opt_names{name_idx});
    end
end

end

