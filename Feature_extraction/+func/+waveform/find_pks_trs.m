function points_out = find_pks_trs(ts,type, plot_flag)
% This function locates the peak or troughs of an input signal
% Code adapted from P.Charlton:
% https://github.com/peterhcharlton/pulse-analyse
%
% INPUT: ts: Signal in
%        type: string indicating whether to locate peaks (pks) or troughs (trs)
%        plot_flag
%
% OUTPUT: points_out: Location of requested points

%pks = peaks. trs = troughs
narginchk(2, inf)
if nargin < 3
    plot_flag = 0;
end

ts_in = ts; % For plotting
if strcmp(type, 'trs') || strcmp(type, 'tr')
    ts = -ts;
elseif ~(strcmp(type, 'pks') || strcmp(type, 'pk'))
    warning('Type not known -- assumming peak')
end


temp1 = find(ts(2:end-1) > ts(1:(end-2)) & ts(2:end-1) > ts(3:end) );

temp2 = find(ts(2:end-2) > ts(1:(end-3)) & ts(2:end-2)== ts(3:(end-1)) & ts(3:end-1) > ts(4:end) );

points_out = unique([temp1; temp2]) +1;
%%
if plot_flag
    figure
    plot(ts_in, 'k')
    hold on
    scatter(points_out, ts_in(points_out), 'ro')
    xlabel('Sample number')
    ylabel('Signal')
    if strcmp(type, 'tr') || strcmp(type, 'trs')
        title('Onset detection')
    else
        title('Peak detection')
    end
    func.plot.tightfig();
end
end
