function plot_indices_and_gauss(PPG, pulse_no, do_gauss_legend, opts)
% make plot of individual beat if needed
narginchk(1, inf)
if nargin < 2 || isempty(pulse_no)
    pulse_no =1;
end
if nargin < 3 || isempty(do_gauss_legend)
    do_gauss_legend = false;
end
if nargin < 4
   opts = struct(); 
end
default_opts.normalise_pw = true;
default_opts.plot_areas = 1;
default_opts.plot_third_deriv =0;
default_opts.plot_pw_only = 0;
default_opts.save_folder =[];
opts = func.aux_functions.update_with_default_opts(default_opts, opts);
%%
warning('off', 'MATLAB:declareGlobalBeforeUse')

params.pulse_no = pulse_no;

%If the indicies dont exist then run get_ppg_indices.m
if ~isfield(PPG, 'fid_pts')
    PPG= func.pulsew.get_ppg_fid_pts(PPG);
end
fid_pts = PPG.fid_pts;
if ~ isfield(fid_pts, 'dic')
    warning('Cannot make plot')
    return
end
%get colours for plotting
global colours

% colours = constants_def('Colours');
colours = func.aux_functions.define_colours();

%% Segment pulse and its derivatives

pulse_quality = PPG.sqi_beat;
pulse_quality = pulse_quality > 0;
if ~pulse_quality(params.pulse_no)
    %find the closest high quality beat
    %1 look at all beats after the current beat
    pulse_quality_temp = pulse_quality;
    pulse_quality_temp( 1:params.pulse_no-1) = 0;
    temp1 = find(pulse_quality_temp & ~isnan(fid_pts.a.ind), 1);
    
    %2 look at all beats before the current beat
    pulse_quality_temp = pulse_quality;
    pulse_quality_temp(params.pulse_no:end) = 0;
    temp2 = find(pulse_quality_temp & ~isnan(fid_pts.a.ind),1, 'last');
    
    [~,idx] = min([abs(temp1 - params.pulse_no), abs(temp2 - params.pulse_no)]);
    if idx ==1
        temp = temp1;
    else
        temp = temp2;
    end
    
    if ~isempty(temp)
        params.pulse_no = temp;
    else
        return
    end
end
curr_els = PPG.onsets(params.pulse_no):PPG.onsets(params.pulse_no+1);

% (either filtered or not)
sigs.v        = PPG.ts(curr_els);
sigs.first    = PPG.derivs.first(curr_els);
sigs.second   = PPG.derivs.second(curr_els);
sigs.third    = PPG.derivs.third(curr_els);

sig_names = fieldnames(sigs);

% - subtract baseline from this pulse wave, and normalise
if opts.normalise_pw
    %subtract baseline
    baseline = linspace(sigs.v(1),sigs.v(end),length(sigs.v));
    sigs.v = sigs.v(:) - baseline(:);
    
    %Normalise
    sigs.v = sigs.v - min(sigs.v);
    scale_factor = max(sigs.v);
    
    for s_idx = 1:length(sig_names)
        sigs.(sig_names{s_idx}) = sigs.(sig_names{s_idx}) / scale_factor;
    end
end

% - Plot fiducial points
sigs.fs = PPG.fs;
%     make_plots(sigs, params.pulse_no, fid_pts, up, 0)

%% - setup
paper_size = 1.125*[400,1050];
paper_size(1) = paper_size(1);
paper_size(2) = paper_size(2);

figure('Position', [2240,64, paper_size])
params.ftsize = 1.5*12; params.lwidth = 2;
sigs.t = [0:length(sigs.v)-1]/sigs.fs;
sigs.t = sigs.t/max(sigs.t);

y_offset = 0.06;
y_inc = 0.187;
n_sub = 5;
x_offset = 0.12;
width = 0.88;
%% - plot sig

h_b(1) = subplot('Position', [x_offset,y_offset+(n_sub-1)*y_inc,width,y_inc-0.01]);
sigs = plot_PPG1(sigs, params, fid_pts, opts);

h_b(2) = subplot('Position', [x_offset,y_offset+(n_sub-2)*y_inc,width,y_inc-0.01]);
sigs = plot_PPG2(sigs, params, fid_pts, opts);


h_b(3) = subplot('Position', [x_offset,y_offset+(n_sub-3)*y_inc,width,y_inc - 0.01]);
sigs = plot_gauss(sigs, params, PPG, do_gauss_legend);


h_b(4) = subplot('Position', [x_offset,y_offset+(n_sub-4)*y_inc,width,y_inc - 0.01]);
sigs = plot_VPG(sigs, params,fid_pts,  opts);

h_b(5) = subplot('Position', [x_offset,y_offset+(n_sub-5)*y_inc,width,y_inc - 0.01]);
plot_APG(sigs, params, fid_pts, opts);
%% - plot Second derivative


linkaxes(h_b, 'x')
end

function normalised  = coords_to_pos(x_coord, y_coord)

pos = get(gca, 'Position');
normalised(1) = (x_coord - min(xlim))/diff(xlim) * pos(3) + pos(1);
normalised(2) = (y_coord - min(ylim))/diff(ylim) * pos(4) + pos(2);

end


function sigs = plot_PPG1(sigs, params, fid_pts, opts)
global colours
% plot baseline curve
plot(sigs.t, sigs.v, 'color', colours.black, 'LineWidth', params.lwidth); hold on,


% plot salient points
% pt_names = {'dia', 'dic', 's', 'halfpoint', 'tangent', 'f1', 'f2'}; %, 'p1in', 'p2in'};
pt_names = {'dia', 'dic', 's', 'f1', 'f2'}; %, 'p1in', 'p2in'};
hspace0 = 0;
for pt_no = 1 : length(pt_names)
    eval(['curr_pt.el = fid_pts.' pt_names{pt_no} '.ind(params.pulse_no) - fid_pts.f1.ind(params.pulse_no)+1;']);
    if isnan(curr_pt.el)
        continue
    end
    
    % plot point
    curr_text = pt_names{pt_no};
    
    
    if any(strcmp({'f1', 'f2'},curr_text))
        eval(['sigs.pts.' curr_text ' = curr_pt.el;']);
        continue
    end
    
    
    
    if strcmp(curr_text, 'tangent')
        curr_pt.v = sigs.v(1);
        curr_pt.t = curr_pt.el/PPG.fs;
    else
        curr_pt.v = sigs.v(curr_pt.el);
        curr_pt.t = sigs.t(curr_pt.el);
    end
    plot(curr_pt.t, curr_pt.v, 'o', 'color', colours.red)
    
    % annotate point
    
    vspace0 = 0.12*range(sigs.v);
    plot_text = curr_text;
    switch curr_text
        case {'dia','p2'}
            hspace0 = 0.04;
            if (strcmp(curr_text, 'p2pk') || strcmp(curr_text, 'p2in')) && curr_pt.el < (fid_pts.s.ind(params.pulse_no) - fid_pts.f1.ind(params.pulse_no)+1)
                hspace0 = -0.04;
            elseif (strcmp(curr_text, 'p1pk') || strcmp(curr_text, 'p1in'))
                hspace0 = 0;
            end
            if strcmp(curr_text, 'dia')
                plot_text = 'D';
            end
        case 'dic'
            hspace0 = -0.01;
            vspace0 = -1*vspace0;
            plot_text = 'N';
        case {'f1', 'p1','p1in'}
            hspace0 = -0.04;
        case 's'
            hspace0 = 0;
            plot_text = 'S';
        case {'halfpoint', 'tangent'}
            plot_text = curr_text(1:4);
            hspace0 = 0;
    end
    
    if strcmp(curr_text, 'tangent')
        text(curr_pt.t+hspace0, curr_pt.v+vspace0 , plot_text,'FontSize', params.ftsize, 'Color', colours.red, 'HorizontalAlignment', 'center');
    else
        text(sigs.t(curr_pt.el)+hspace0, sigs.v(curr_pt.el)+vspace0 , plot_text,'FontSize', params.ftsize, 'Color', colours.red, 'HorizontalAlignment', 'center');
    end
    eval(['sigs.pts.' curr_text ' = curr_pt.el;']);
end

% set limits
curr_range = range(sigs.v);
% if plot_inds
factor1 = 0.4; factor2 = 0.2;
% else
%     factor1 = 0.15; factor2 = 0.15;
% end
% factor1 = 1;
ylims = [min(sigs.v)-factor2*curr_range, max(sigs.v)+factor1*curr_range];

ylim(ylims)
curr_range = range(sigs.t);
xlims = [sigs.t(1)-0.05*curr_range, sigs.t(end)+0.1*curr_range];
xlim_labs = 0:0.25:0.25*ceil(sigs.t(end)*4);
xlim(xlims)

% set labels
ylab = ylabel('PPG', 'FontSize', params.ftsize, 'Rotation', 90);
%     set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(gca, 'FontSize', params.ftsize -2, 'XTick', xlim_labs, 'XGrid', 'on')
if ~opts.normalise_pw
    ytick_vals = [round(min(sigs.v),3,'significant'), round(max(sigs.v),3,'significant')];
else
    ytick_vals = [];
end
set(gca, 'YTick', ytick_vals)
box off
if opts.plot_pw_only
    xlabel('Time (s)', 'FontSize', params.ftsize)
else
    set(gca, 'XTickLabel', {})
end

% Plot indices

color = 0.4;

% - delta T
ht_annot = ylims(2)-0.11*range(ylims);
plot(sigs.t(sigs.pts.s)*[1,1], [ht_annot, sigs.v(sigs.pts.s)], '--', 'color', color*ones(1,3))
plot(sigs.t(sigs.pts.dia)*[1,1], [ht_annot, sigs.v(sigs.pts.dia)], '--', 'color', color*ones(1,3))
normalised1  = coords_to_pos(sigs.t(sigs.pts.s), ht_annot);
normalised2  = coords_to_pos(sigs.t(sigs.pts.dia), ht_annot);
ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
ah.Head1Length = 7; ah.Head2Length = 7;
ah.Head1Width = 7; ah.Head2Width = 7;
new_ht_annot= ht_annot +0.01*range(ylims);
text(mean([sigs.t(sigs.pts.s),sigs.t(sigs.pts.dia)]), new_ht_annot, '\DeltaT','FontSize', params.ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% - Crest Time
ht_annot = ylims(2)-0.11*range(ylims);
plot(sigs.t(1)*[1,1], [ht_annot, sigs.v(1)], '--', 'color', color*ones(1,3))
plot(sigs.t(sigs.pts.s)*[1,1], [ht_annot, sigs.v(sigs.pts.s)], '--', 'color', color*ones(1,3))
normalised1  = coords_to_pos(sigs.t(1), ht_annot);
normalised2  = coords_to_pos(sigs.t(sigs.pts.s), ht_annot);
ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
ah.Head1Length = 7; ah.Head2Length = 7;
ah.Head1Width = 7; ah.Head2Width = 7;
new_ht_annot= ht_annot +0.01*range(ylims);
text(mean([sigs.t(1),sigs.t(sigs.pts.s)]), new_ht_annot, 'CT','FontSize', params.ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% - Reflection Index
t_annot = sigs.t(end)+0.05;
plot([sigs.t(end), t_annot], [1,1]*sigs.v(end), '--', 'color', color*ones(1,3))
plot([sigs.t(sigs.pts.dia), t_annot], [1,1]*sigs.v(sigs.pts.dia), '--', 'color', color*ones(1,3))
normalised1  = coords_to_pos(t_annot, sigs.v(end));
normalised2  = coords_to_pos(t_annot, sigs.v(sigs.pts.dia));
ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
ah.Head1Length = 7; ah.Head2Length = 7;
ah.Head1Width = 7; ah.Head2Width = 7;
new_ht_annot= mean([sigs.v(sigs.pts.dia), sigs.v(end)]);
text(sigs.t(end), new_ht_annot, 'RI','FontSize', params.ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% - Systolic time
ht_annot = ylims(1)+0.02*range(ylims);
plot(sigs.t(1)*[1,1], [ht_annot, sigs.v(1)], '--', 'color', color*ones(1,3))
plot(sigs.t(sigs.pts.dic)*[1,1], [ht_annot, sigs.v(sigs.pts.dic)], '--', 'color', color*ones(1,3))
normalised1  = coords_to_pos(sigs.t(1), ht_annot);
normalised2  = coords_to_pos(sigs.t(sigs.pts.dic), ht_annot);
ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
ah.Head1Length = 7; ah.Head2Length = 7;
ah.Head1Width = 7; ah.Head2Width = 7;
new_ht_annot= ht_annot +0.01*range(ylims);
% text(mean([sigs.t(1),sigs.t(sigs.pts.dic)]), new_ht_annot, 'Systole','FontSize', params.ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(mean([sigs.t(1),sigs.t(sigs.pts.dic)]), new_ht_annot, 'T_{Sys}','FontSize', params.ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% - Diastolic time
ht_annot = ylims(1)+0.02*range(ylims);
plot(sigs.t(end)*[1,1], [ht_annot, sigs.v(1)], '--', 'color', color*ones(1,3))
plot(sigs.t(sigs.pts.dic)*[1,1], [ht_annot, sigs.v(sigs.pts.dic)], '--', 'color', color*ones(1,3))
normalised1  = coords_to_pos(sigs.t(sigs.pts.dic), ht_annot);
normalised2  = coords_to_pos(sigs.t(end), ht_annot);
ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
ah.Head1Length = 7; ah.Head2Length = 7;
ah.Head1Width = 7; ah.Head2Width = 7;
new_ht_annot= ht_annot +0.01*range(ylims);
% text(mean([sigs.t(sigs.pts.dic),sigs.t(end)]), new_ht_annot, 'Diastole','FontSize', params.ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(mean([sigs.t(sigs.pts.dic),sigs.t(end)]), new_ht_annot, 'T_{Dias}','FontSize', params.ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');




end


function sigs = plot_PPG2(sigs, params, fid_pts, opts)
global colours
% plot baseline curve
plot(sigs.t, sigs.v, 'color', colours.black, 'LineWidth', params.lwidth); hold on,
% plot(sigs.t, sigs.v, 'LineWidth', params.lwidth); hold on,
if  opts.plot_areas
    
    % Areas: Systolic
    el = fid_pts.dic.ind(params.pulse_no) - fid_pts.f1.ind(params.pulse_no)+1;
    rel_els = 1:el-1;
    offset = 0.05*range(sigs.v);
    %     h = fill([sigs.t(rel_els), sigs.t(rel_els(end))], [sigs.v(rel_els); sigs.v(rel_els(1))], [1,0.8,0.8]);
    h = fill([sigs.t(rel_els), sigs.t(rel_els(end))], [sigs.v(rel_els); sigs.v(rel_els(1))], colours.yellow);
    h.LineStyle = 'none';
    h.FaceAlpha = 0.5;
    
    % Areas: Diastolic
    el2 = fid_pts.f2.ind(params.pulse_no) - fid_pts.f1.ind(params.pulse_no)+1;
    rel_els = el+1:el2;
    offset = 0.05*range(sigs.v);
    %     h = fill([sigs.t(rel_els), sigs.t(rel_els(1))], [sigs.v(rel_els); sigs.v(rel_els(end))], [0.8,0.8,1]);
    h = fill([sigs.t(rel_els), sigs.t(rel_els(1))], [sigs.v(rel_els); sigs.v(rel_els(end))], colours.green);
    h.LineStyle = 'none';
    h.FaceAlpha = 0.5;
    plot(sigs.t, sigs.v,'color', colours.black, 'LineWidth', params.lwidth)
    
end


% plot salient points
% pt_names = { 'dic', 'f1', 'f2'}; %, 'p1in', 'p2in'};
pt_names = {'dic'}; %, 'p1in', 'p2in'};
hspace0 = 0;
for pt_no = 1 : length(pt_names)
    eval(['curr_pt.el = fid_pts.' pt_names{pt_no} '.ind(params.pulse_no) - fid_pts.f1.ind(params.pulse_no)+1;']);
    if isnan(curr_pt.el)
        continue
    end
    
    % plot point
    curr_text = pt_names{pt_no};
    
    if any(strcmp({'f1', 'f2'},curr_text))
        eval(['sigs.pts.' curr_text ' = curr_pt.el;']);
        continue
    end
    
    
    
    if strcmp(curr_text, 'tangent')
        curr_pt.v = sigs.v(1);
        curr_pt.t = curr_pt.el/PPG.fs;
    else
        curr_pt.v = sigs.v(curr_pt.el);
        curr_pt.t = sigs.t(curr_pt.el);
    end
    plot(curr_pt.t, curr_pt.v, 'o', 'color', colours.red)
    
    % annotate point
    
    vspace0 = 0.12*range(sigs.v);
    plot_text = curr_text;
    switch curr_text
        case {'dia','p2'}
            hspace0 = 0.04;
            if (strcmp(curr_text, 'p2pk') || strcmp(curr_text, 'p2in')) && curr_pt.el < (fid_pts.s.ind(params.pulse_no) - fid_pts.f1.ind(params.pulse_no)+1)
                hspace0 = -0.04;
            elseif (strcmp(curr_text, 'p1pk') || strcmp(curr_text, 'p1in'))
                hspace0 = 0;
            end
        case 'dic'
            hspace0 = -0.01;
            vspace0 = 1*vspace0;
            plot_text ='N';
        case {'f1', 'p1','p1in'}
            hspace0 = -0.04;
        case 's'
            hspace0 = 0;
        case {'halfpoint', 'tangent'}
            plot_text = curr_text(1:4);
            hspace0 = 0;
    end
    
    if strcmp(curr_text, 'tangent')
        text(curr_pt.t+hspace0, curr_pt.v+vspace0 , plot_text,'FontSize', params.ftsize, 'Color', colours.red, 'HorizontalAlignment', 'center');
    else
        text(sigs.t(curr_pt.el)+hspace0, sigs.v(curr_pt.el)+vspace0 , plot_text,'FontSize', params.ftsize, 'Color', colours.red, 'HorizontalAlignment', 'center');
    end
    eval(['sigs.pts.' curr_text ' = curr_pt.el;']);
end

% set limits
curr_range = range(sigs.v);
% if plot_inds
factor1 = 0.4; factor2 = 0.2;
% else
%     factor1 = 0.15; factor2 = 0.15;
% end
% factor1 = 1;
ylims = [min(sigs.v)-factor2*curr_range, max(sigs.v)+factor1*curr_range];

ylim(ylims)
curr_range = range(sigs.t);
xlims = [sigs.t(1)-0.05*curr_range, sigs.t(end)+0.1*curr_range];
xlim_labs = 0:0.25:0.25*ceil(sigs.t(end)*4);
xlim(xlims)

% set labels
ylab = ylabel('PPG', 'FontSize',params.ftsize, 'Rotation', 90);
%     set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(gca, 'FontSize', params.ftsize -2, 'XTick', xlim_labs, 'XGrid', 'on')
if ~opts.normalise_pw
    ytick_vals = [round(min(sigs.v),3,'significant'), round(max(sigs.v),3,'significant')];
else
    ytick_vals = [];
end
set(gca, 'YTick', ytick_vals)
box off
if opts.plot_pw_only
    xlabel('Normalised time', 'FontSize', params.ftsize)
else
    set(gca, 'XTickLabel', {})
end




% - Systolic area
% ht_annot = mean([sigs.v(1), sigs.v(sigs.pts.s)]);
ht_annot = sigs.v(sigs.pts.f1) + 0.1;
text(mean([sigs.t(sigs.pts.dic),sigs.t(1)]), ht_annot, 'A1','FontSize', params.ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

% - Diastolic area
% ht_annot = mean([sigs.v(1), sigs.v(sigs.pts.s)]);
text(mean([sigs.t(sigs.pts.dic),sigs.t(end)]), ht_annot, 'A2','FontSize', params.ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');


%STT
col = colours.purple;
plot(sigs.t([sigs.pts.f1,sigs.pts.s]), sigs.v([sigs.pts.f1, sigs.pts.s]), 'color', col, 'LineWidth', params.lwidth),
curr_range = range(sigs.t);
x = mean(sigs.t([sigs.pts.f1,sigs.pts.s]));
y = mean(sigs.v([sigs.pts.f1,sigs.pts.f1,sigs.pts.s,sigs.pts.s]));
text(x-0.03*curr_range, y, 'STT','FontSize', params.ftsize, 'Color', col, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Rotation', 50);



find_crossing = @(v, t) find((v(:)-t).*circshift((v(:)-t), [-1 0]) <= 0);
% Width 25
cutoff_25 = 0.25*(sigs.v(sigs.pts.s) - sigs.v(sigs.pts.f1));
cutoff_25 = cutoff_25 + sigs.v(sigs.pts.f1);
crossing_25 = find_crossing(sigs.v, cutoff_25);
% plot(sigs.t([crossing_25(1),crossing_25(end)]),[cutoff_25 cutoff_25], '--', 'color', color*ones(1,3))


% ht_annot = ylims(1)+0.02*range(ylims);
params.ftsize_width = 12;
ht_annot = cutoff_25;
normalised1  = coords_to_pos(sigs.t(crossing_25(1)), ht_annot);
normalised2  = coords_to_pos(sigs.t(crossing_25(end)), ht_annot);
ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)], 'Color',  0.2*ones(1,3));
ah.Head1Length = 7; ah.Head2Length = 7;
ah.Head1Width = 7; ah.Head2Width = 7;
new_ht_annot= ht_annot -0.1*range(ylims);
text(sigs.t(sigs.pts.dic), new_ht_annot, 'Width 25','FontSize', params.ftsize_width, 'Color',  0.2*ones(1,3), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');




% Width 50
cutoff_50 = 0.5*(sigs.v(sigs.pts.s) - sigs.v(sigs.pts.f1));
cutoff_50 = cutoff_50 + sigs.v(sigs.pts.f1);
crossing_50 = find_crossing(sigs.v, cutoff_50);
% plot(sigs.t([crossing_50(1),crossing_50(end)]),[cutoff_50 cutoff_50], '--', 'color', color*ones(1,3))
ht_annot = cutoff_50;
normalised1  = coords_to_pos(sigs.t(crossing_50(1)), ht_annot);
normalised2  = coords_to_pos(sigs.t(crossing_50(end)), ht_annot);
ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)], 'Color',  0.2*ones(1,3));
ah.Head1Length = 7; ah.Head2Length = 7;
ah.Head1Width = 7; ah.Head2Width = 7;
new_ht_annot= ht_annot -0.1*range(ylims);
text(sigs.t(sigs.pts.dic), new_ht_annot, 'Width 50','FontSize', params.ftsize_width, 'Color', 0.2*ones(1,3), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

end

function sigs = plot_gauss(sigs, params, PPG,do_legend)
global colours
hppg = plot(sigs.t, sigs.v, 'color', 'k', 'LineWidth', params.lwidth); hold on,
gaussian = @(b,x) b(1) * exp(-(x - b(2)).^2/b(3));
modelfun = @(b,x) gaussian(b(1:3), x(:,1)) + gaussian(b(4:6), x(:,1))+ gaussian(b(7:9), x(:,1))+ gaussian(b(10:12), x(:,1));


% params = table2array(gaussian_parameters);
gauss_params = [PPG.norm_fid_pts.g1.amp(params.pulse_no);
    PPG.norm_fid_pts.g1.mu(params.pulse_no);
    PPG.norm_fid_pts.g1.sigma(params.pulse_no);
    PPG.norm_fid_pts.g2.amp(params.pulse_no);
    PPG.norm_fid_pts.g2.mu(params.pulse_no);
    PPG.norm_fid_pts.g2.sigma(params.pulse_no);
    PPG.norm_fid_pts.g3.amp(params.pulse_no);
    PPG.norm_fid_pts.g3.mu(params.pulse_no);
    PPG.norm_fid_pts.g3.sigma(params.pulse_no);
    PPG.norm_fid_pts.g4.amp(params.pulse_no);
    PPG.norm_fid_pts.g4.mu(params.pulse_no);
    PPG.norm_fid_pts.g4.sigma(params.pulse_no);
    ];

tot = modelfun(gauss_params, sigs.t');
g1 = gaussian(gauss_params(1:3), sigs.t');
g2 = gaussian(gauss_params(4:6), sigs.t');
g3 = gaussian(gauss_params(7:9), sigs.t');
g4 = gaussian(gauss_params(10:12), sigs.t');
sys = g1+ g2;

%     h = fill([0, sigs.t], [0; sys], colours.green);
%     h.LineStyle = 'none';
%     h.FaceAlpha = 0.2;
%     % hatchfill(h, 'single',45,5,'none');
spacing = 10;
%     hatchfill(h, 'single',45,spacing,colours.green);


%     h = fill([0, sigs.t], [0; g3], colours.yellow);
%     h.LineStyle = 'none';
%     h.FaceAlpha = 0.2;
%     hatchfill(h, 'single',45,spacing,colours.yellow);


h1 = plot(sigs.t, tot, '--', 'Color', colours.blue, 'LineWidth', params.lwidth);
h2 = plot(sigs.t,sys , ':', 'Color', colours.green, 'LineWidth', params.lwidth);


h3 = plot(sigs.t, g1, '--', 'Color', colours.bright_blue, 'LineWidth', params.lwidth/2);
h4 = plot(sigs.t, g2, '--', 'Color', colours.dark_blue, 'LineWidth', params.lwidth/2);
h5 = plot(sigs.t, g3, '--', 'Color', colours.yellow, 'LineWidth', params.lwidth/2);
h6 = plot(sigs.t, g4, '--', 'Color', colours.red, 'LineWidth', params.lwidth/2);

%     labels = {'PPG', 'Model', 'Sys', 'g1', 'g2', 'g3', 'g4'};
labels = {'PPG', 'Model', 'g1', 'g2', 'g3', 'g4'};

set(gca,'xticklabel',{[]})
box off

%% Do annotations
color = 0.15;

factor2 = 0.15; factor1 = 0.1;
curr_range= range(sigs.v);
ylims = [min(sigs.v)-factor2*curr_range, max(sigs.v)+factor1*curr_range];


ht_annot = ylims(2)-0.02*range(ylims);
plot(gauss_params(2)*[1,1], [ht_annot, gauss_params(1)], '--', 'color', color*ones(1,3))
plot(gauss_params(8)*[1,1], [ht_annot, gauss_params(7)], '--', 'color', color*ones(1,3))
normalised1  = coords_to_pos(gauss_params(2), ht_annot);
normalised2  = coords_to_pos(gauss_params(8), ht_annot);
ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
ah.Head1Length = 7; ah.Head2Length = 7;
ah.Head1Width = 7; ah.Head2Width = 7;
new_ht_annot= ht_annot ;
text((gauss_params(2) + gauss_params(8))/2, new_ht_annot, 'Gauss RTT','FontSize', params.ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');



t = sigs.t(end);
[~,i] = max(g1);
t1 = sigs.t(i);
% [~,i] = max(g2);
% t2 = sigs.t(i);
[amp1,i] = max(sys);
t2 = sigs.t(i);
% amp1 = gauss_params(1);
amp2 = gauss_params(4);

t_annot = t-0.05;
plot([t1, t_annot], [1,1]*amp1, '--', 'color', color*ones(1,3))
plot([t2, t_annot], [1,1]*amp2, '--', 'color', color*ones(1,3))
normalised1  = coords_to_pos(t_annot, amp1);
normalised2  = coords_to_pos(t_annot, amp2);
ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
ah.Head1Length = 7; ah.Head2Length = 7;
ah.Head1Width = 7; ah.Head2Width = 7;
new_ht_annot= (amp1 + amp2)/2;
text(sigs.t(end), new_ht_annot, {'Gauss', 'AI'},'FontSize', params.ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

ylabel('Normalised amplitude')

% set(gca, 'XTickLabel', {})
% set(gca, 'YTickLabel', {})
% xticks = get(gca,'xtick');
% newlabels = arrayfun(@(x) sprintf('%.1f', x./time_factor ), xticks, 'un', 0);
% set(gca,'xticklabel',newlabels);

% legend([hppg, h1, h2, h3, h4, h5, h6], labels, 'Location', 'bestoutside')

ylab = ylabel({'PPG Gauss'}, 'FontSize', params.ftsize, 'Rotation', 90);
%     set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
xlims = [sigs.t(1)-0.05*curr_range, sigs.t(end)+0.1*curr_range];
xlim_labs = 0:0.25:0.25*ceil(sigs.t(end)*4);

set(gca, 'FontSize', params.ftsize -2, 'YTick', [], 'XTick', xlim_labs, 'XGrid', 'on')
%%
set(gca, 'YTick', []);
ylim([ylims(1), 1.3])
%     ylim([0, 1.3])
if do_legend
    legend([hppg, h1, h3, h4, h5, h6], labels)
end
end


function sigs = plot_VPG(sigs, params,fid_pts,  opts)
global colours
% Plot x-axis
plot([-10, 10], [0,0], '--k'); hold on

% plot baseline curve
plot(sigs.t, sigs.first, 'color', colours.black, 'LineWidth', params.lwidth); hold on,

% plot salient points
pt_names = {'W', 'dic'};
curr_range = range(sigs.first);
vspace0 = 0.1*curr_range;
hspace0 = 0.05;
for pt_no = 1 : length(pt_names)
    eval(['curr_fid_pt = fid_pts.' pt_names{pt_no} '.ind(params.pulse_no);'])
    if isnan(curr_fid_pt)
        continue
    end
    curr_pt.el = curr_fid_pt - fid_pts.f1.ind(params.pulse_no)+1;
    if isnan(curr_pt.el)
        continue
    end
    
    % plot point
    curr_pt.v = sigs.first(curr_pt.el);
    curr_pt.t = sigs.t(curr_pt.el);
    plot(curr_pt.t, curr_pt.v, 'o', 'color', colours.red)
    
    % annotate point
    curr_text = pt_names{pt_no};
    plot_text = curr_text;
    if strcmp(plot_text, 'dic')
        plot_text ='N_+_1';
    end
    text(sigs.t(curr_pt.el)+hspace0, sigs.first(curr_pt.el)+vspace0 , plot_text,'FontSize', params.ftsize, 'Color', colours.red, 'HorizontalAlignment', 'center');
    eval(['sigs.pts.' curr_text ' = curr_pt.el;']);
    
end

% set limits
curr_range = range(sigs.t);
xlims = [sigs.t(1)-0.05*curr_range, sigs.t(end)+0.1*curr_range];
xlim_labs = 0:0.25:0.25*ceil(sigs.t(end)*4);
xlim(xlims)
curr_range = range(sigs.first);
ylim([min(sigs.first)-0.05*curr_range, max(sigs.first)+0.15*curr_range]); ylims = ylim;

% set labels
ylab = ylabel({'VPG'}, 'FontSize', params.ftsize, 'Rotation', 90);
% set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(gca, 'FontSize', params.ftsize -2, 'YTick', [])
set(gca, 'XTick', xlim_labs, 'XTickLabel', {}, 'XGrid', 'on')
if ~opts.normalise_pw
    ytick_vals = [0, round(max(sigs.first),3,'significant')];
    set(gca, 'YTick', ytick_vals)
else
    set(gca, 'YTick', 0);
end
box off

% Plot indices

% - ms
ht_annot = mean([sigs.first(sigs.pts.W), min(sigs.first)]);
[~, temp] = min(sigs.first(sigs.pts.W:sigs.pts.dic));
el = temp + sigs.pts.W -1;
%         plot([sigs.t(sigs.pts.W), sigs.t(el)], sigs.first(sigs.pts.W)*[1,1], '--', 'color', color*ones(1,3))
%     plot([sigs.t(sigs.pts.ms), sigs.t(el)], sigs.first(el)*[1,1], '--', 'color', color*ones(1,3))
%         normalised1  = coords_to_pos(sigs.t(el), 0);
%         normalised2  = coords_to_pos(sigs.t(el), sigs.first(sigs.pts.W));
%         ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
%         curr_range = range(sigs.t);
%         text(sigs.t(el)+0.06*curr_range, ht_annot, 'W','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

end

function plot_APG(sigs, params, fid_pts, opts)
global colours

% Plot x-axis
plot([-10, 10], [0,0], '--k'); hold on

% plot baseline curve
curr_color = 0.4*0.5*[1,2,1];
plot(sigs.t, sigs.second, 'color', colours.black, 'LineWidth', params.lwidth); hold on,

% plot salient points
pt_names = {'a', 'b', 'c', 'd', 'e', 'f'};
curr_range = range(sigs.second);
vspace_const = 0.08;
hspace0 = 0;
for pt_no = 1 : length(pt_names)
    eval(['curr_pt.el = fid_pts.' pt_names{pt_no} '.ind(params.pulse_no) - fid_pts.f1.ind(params.pulse_no)+1;']);
    if isnan(curr_pt.el)
        continue
    end
    
    % plot point
    curr_pt.v = sigs.second(curr_pt.el);
    curr_pt.t = sigs.t(curr_pt.el);
    plot(curr_pt.t, curr_pt.v, 'o', 'color', colours.red)
    
    % annotate point
    curr_text = pt_names{pt_no};
    switch curr_text
        case {'a','c', 'e'}
            vspace0 = (1.3*vspace_const)*curr_range;
        case {'b', 'd', 'f'}
            vspace0 = -1.3*vspace_const*curr_range;
    end
    text(sigs.t(curr_pt.el), sigs.second(curr_pt.el)+vspace0 , curr_text,'FontSize', params.ftsize, 'Color', colours.red, 'HorizontalAlignment', 'center');
    eval(['sigs.pts.' curr_text ' = curr_pt.el;']);
end

% set limits
curr_range = range(sigs.t);
xlims = [sigs.t(1)-0.05*curr_range, sigs.t(end)+0.1*curr_range];
xlim_labs = 0:0.25:0.25*ceil(sigs.t(end)*4);
xlim(xlims)
xlim(xlims)
curr_range = range(sigs.second);
ylim([min(sigs.second)-0.15*curr_range, max(sigs.second)+0.15*curr_range])

% set labels
ylab = ylabel({'APG'}, 'FontSize', params.ftsize, 'Rotation', 90);
%     set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(gca, 'FontSize', params.ftsize -2, 'YTick', [], 'XTick', xlim_labs, 'XGrid', 'on')
box off
if ~opts.plot_third_deriv
    xlabel('Normalised time', 'FontSize', params.ftsize)
else
    set(gca, 'XTickLabel', {})
end
if ~opts.normalise_pw
    ytick_vals = [round(min(sigs.second),3,'significant'), 0, round(max(sigs.second),3,'significant')];
    set(gca, 'YTick', ytick_vals)
else
    set(gca, 'YTick', 0);
end



% - slope b-c
plot(sigs.t([sigs.pts.b,sigs.pts.c]), sigs.second([sigs.pts.b, sigs.pts.c]), 'k', 'LineWidth', params.lwidth),
curr_range = range(sigs.t);
x = mean(sigs.t([sigs.pts.b,sigs.pts.c]));
y = mean(sigs.second([sigs.pts.b,sigs.pts.b,sigs.pts.b,sigs.pts.c]));
text(x-0.1*curr_range, y, 'slope_{b-c}','FontSize', params.ftsize-3, 'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');


% - slope b-d
plot(sigs.t([sigs.pts.b,sigs.pts.d]), sigs.second([sigs.pts.b, sigs.pts.d]), 'k', 'LineWidth', params.lwidth),
curr_range = range(sigs.t);
x = mean(sigs.t([sigs.pts.b,sigs.pts.d]));
y = mean(sigs.second([sigs.pts.b,sigs.pts.b,sigs.pts.b,sigs.pts.d]));
text(x+0.03*curr_range, y, 'slope_{b-d}','FontSize', params.ftsize-3, 'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end