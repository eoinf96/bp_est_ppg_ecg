function opts = update_with_default_opts(opts,default_opts)

opt_names = fieldnames(default_opts); 
for name_idx = 1:length(opt_names)
    if ~isfield(opts,opt_names{name_idx} )
        opts.(opt_names{name_idx}) = default_opts.(opt_names{name_idx});
    end
end

end

