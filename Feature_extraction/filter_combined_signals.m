function [data_out ] = filter_combined_signals(sig_in, filter_config, remove_powerline)
    narginchk(2, inf)
    if nargin < 3
       remove_powerline = false; 
    end
    %%   
    
    data_out = sig_in.ts;
    if isempty(data_out)
       return 
    end
    
    if ~any(isnan(sig_in.ts))
        if remove_powerline
            sig_in.ts =  pe.sigproc.ecg.remove_powerline_noise(sig_in.ts, sig_in.fs,[], filter_config.mains_f);
        end
        
        data_out =  pe.sigproc.filter.filter(sig_in.ts, sig_in.fs, filter_config);
        return
    end
    exit = 0;
    good_locs = double(~isnan(sig_in.ts));
    while ~exit
        starting_loc = find(good_locs, 1, 'first');
        if isempty(starting_loc)
           exit = 1;
           continue
        end
        good_locs(1:starting_loc) = 1;
        end_loc = find(1-good_locs, 1, 'first');
        if isempty(end_loc)
            end_loc = length(good_locs)+1;
            exit = 1;
        end
        end_loc = end_loc -1;
        good_locs(1:end_loc) = 0;
        if remove_powerline
            sig_in.ts(starting_loc:end_loc) =  pe.sigproc.ecg.remove_powerline_noise(sig_in.ts(starting_loc:end_loc), sig_in.fs,[], filter_config.mains_f);
        end
        data_out(starting_loc:end_loc) =  pe.sigproc.filter.filter(sig_in.ts(starting_loc:end_loc), sig_in.fs, filter_config);        
    end
    
end