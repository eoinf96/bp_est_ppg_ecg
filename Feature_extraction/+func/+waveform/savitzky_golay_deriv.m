function deriv = savitzky_golay_deriv(ts, deriv_order, win_size)
% This function computes the Savitzky-Golay derivative filter using filter
% coefficients
%
% INPUT: ts: Signal in
%        deriv_order: Order of derivative
%        win_size: Window size
%
% OUTPUT: deriv: Derivative out
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


switch deriv_order
    case 0
        % - smoothing
        switch win_size
            case 5
                coeffs = [-3, 12, 17, 12, -3];
                norm_factor = 35;
            case 7
                coeffs = [-2, 3, 6, 7, 6, 3, -2];
                norm_factor = 21;
            case 9
                coeffs = [-21, 14, 39, 54, 59, 54, 39, 14, -21];
                norm_factor = 231;
            otherwise
                error('Can''t do this window size')
        end
    case 1
        % - first derivative
        switch win_size
            case 5
                coeffs = -2:2;
                norm_factor = 10;
            case 7
                coeffs = -3:3;
                norm_factor = 28;
            case 9
                coeffs = -4:4;
                norm_factor = 60;
            otherwise
                error('Can''t do this window size')
        end
        
    case 2
        % - second derivative
        switch win_size
            case 5
                coeffs = [2,-1,-2,-1,2];
                norm_factor = 7;
            case 7
                coeffs = [5,0,-3,-4,-3,0,5];
                norm_factor = 42;
            case 9
                coeffs = [28,7,-8,-17,-20,-17,-8,7,28];
                norm_factor = 462;
            otherwise
                error('Can''t do this window size')
        end
        
    case 3
        % - third derivative
        switch win_size
            case 5
                coeffs = [-1,2,0,-2,1];
                norm_factor = 2;
            case 7
                coeffs = [-1,1,1,0,-1,-1,1];
                norm_factor = 6;
            case 9
                coeffs = [-14,7,13,9,0,-9,-13,-7,14];
                norm_factor = 198;
            otherwise
                error('Can''t do this window size')
        end
        
    case 4
        % - fourth derivative
        switch win_size
            case 7
                coeffs = [3,-7,1,6,1,-7,3];
                norm_factor = 11;
            case 9
                coeffs = [14,-21,-11,9,18,9,-11,-21,14];
                norm_factor = 143;
            otherwise
                error('Can''t do this window size')
        end
        
    otherwise
        error('Can''t do this order of derivative')
end

if rem(deriv_order, 2) == 1
    coeffs = -1*coeffs;
end

A = [1,0];
filtered_sig = filter(coeffs, A, ts);
s=length(ts);
half_win_size = floor(win_size*0.5);
deriv=[filtered_sig(win_size)*ones(half_win_size,1);filtered_sig(win_size:s);filtered_sig(s)*ones(half_win_size,1)];
deriv = deriv/norm_factor;

end