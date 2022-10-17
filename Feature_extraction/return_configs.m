function configs = return_configs()

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

