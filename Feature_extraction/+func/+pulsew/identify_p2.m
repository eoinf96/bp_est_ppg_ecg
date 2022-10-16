function temp_p2 = identify_p2(curr, temp_d, temp_p1, temp_dic)

% Find "p2" from the minimum value of the third derivative immediately before "d"
td_trs = func.waveform.find_pks_trs(curr.derivs.third, 'tr');
temp = find(td_trs < temp_d,1,'last'); clear d_pks
temp_p2 = td_trs(temp);

% unless c=d, in which case p2 will now be before p1, so take the minimum
% value of the third derivative immediately after "d"
if temp_p2 < temp_p1
    temp_p2 = td_trs(find(td_trs<temp_dic,1,'last'));
end

% check to see if there is a maximum in the signal between the current
% estimate of p2, and temp_dic. If so, take this instead
pks = func.waveform.find_pks_trs(curr.ts, 'pk');
try
temp = find(pks> temp_p2 & pks < temp_dic);
catch
    temp_p2 = nan;
    return
end

if length(temp) == 1
    temp_p2 = pks(temp);
elseif length(temp) == 2
    temp_p2 = pks(temp(2));
elseif length(temp) > 1
%     fprintf('\nCheck this - identify p2')
end

end