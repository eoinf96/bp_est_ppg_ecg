function trs = find_pks_trs(sig,type, plot_flag)
%pks = peaks. trs = troughs 
narginchk(2, inf)
if nargin < 3
    plot_flag = 0;
end

if strcmp(type, 'trs') || strcmp(type, 'tr')
    sig = -sig;
elseif ~(strcmp(type, 'pks') || strcmp(type, 'pk'))
    warning('Type not known -- assumming peak')
end


temp1 = find(sig(2:end-1) > sig(1:(end-2)) & sig(2:end-1) > sig(3:end) );

temp2 = find(sig(2:end-2) > sig(1:(end-3)) & sig(2:end-2)== sig(3:(end-1)) & sig(3:end-1) > sig(4:end) );

temp = unique([temp1; temp2]);

trs = temp+1;

%%
if plot_flag
    %get_colours for plotting
    constants_def;
    
    figure
    plot(sig)
    hold on
    scatter(trs, sig(trs))
    xlabel('Sample number')
    ylabel('Signal')
    if strcmp(type, 'tr')
        title('Onset detection')
    else
        title('Peak detection')
    end
    
end
end


%old method for finding peaks and troughs from adjust_fid_points.m

%%%% dsig2_dt2.ts = diff(dsig_dt.ts);
% dsig2_dt2.t = Signal.t(3:end); %Backards difference derivative
% dsig2_dt2.ts = filtfilt(b,a,dsig2_dt2.ts); %Zero phase filtering
% 
% %Zero crossing function:: this function locates the index to the left of where we go from negative to positive 
% find_zero_crossing = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
% zc = find_zero_crossing(dsig_dt.ts); %Find zero crossings of derivative for max and minimum detection
% zc(diff(zc) == 1) = [];
% 
% %As we cannot know the second derivative of the last point 
% zc(zc ==length(dsig_dt.t)) = [];
% 
% idx_onsets =1;
% idx_peaks = 1;
% %preallocate for memory
% zc_onsets  = nan(size(zc)); 
% zc_peaks  = nan(size(zc));
% for idx = 1:length(zc)
%    if dsig2_dt2.ts(zc(idx)) > 0 
%         %point is onset
%         zc_onsets(idx_onsets) = zc(idx) + 1; % as dppg/dt is backwards difference
%         idx_onsets = idx_onsets +1;
%    elseif dsig2_dt2.ts(zc(idx)) < 0 
%         %point is peak
%         zc_peaks(idx_peaks) = zc(idx) +1;
%         idx_peaks = idx_peaks +1;
% 
%    end
%    % else is stationary point and therefore not a peak or an onset
% end
% zc_onsets(isnan(zc_onsets)) = [];
% zc_peaks(isnan(zc_peaks)) = [];
% 
