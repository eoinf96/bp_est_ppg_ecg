function temp_c = identify_c(curr, temp_b, temp_e)

% Identify C as the highest local maximum on the second derivative between "b" and "e"
pks = func.waveform.find_pks_trs(curr.derivs.second, 'pk');
temp = find(pks > temp_b & pks < temp_e);
[~, rel_tr_el] = max(curr.derivs.second(pks(temp)));
temp_c = pks(temp(rel_tr_el)); clear temp rel_tr_el pks

% If there aren't any peaks that satisfy this criterion ...
if isempty(temp_c)
    % then identify C as the lowest local minimum on the third derivative
    % after "b" and before "e"
    trs = func.waveform.find_pks_trs(curr.derivs.third, 'tr');
    temp = find(trs > temp_b & trs < temp_e);
    [~, rel_tr_el] = min(curr.derivs.third(trs(temp)));
    if ~isempty(rel_tr_el)
        temp_c = trs(temp(rel_tr_el)); clear temp rel_tr_el trs
    end
end

end