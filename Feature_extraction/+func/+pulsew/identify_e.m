function temp_e = identify_e(curr, temp_ms, temp_b)
    

% Find local maxima in the second derivative
pks = func.waveform.find_pks_trs(curr.derivs.second, 'pk');
% Set an upper bound of 60 % of the PW duration
upper_bound = 0.6*length(curr.ts);   % const from: https://en.wikipedia.org/wiki/QT_interval#/media/File:QT_interval_corrected_for_heart_rate.png
%%%%%
% QT interval represents time from ventricular contraction to relaxation.
% Time to dicrotic notch can also be thought of in a similar manner and so
% the two are related
% 0.6 is essentially an upper estimate of the possible duration of systole
% - In MOLLIE for example (and in the figure referenced above), t sys is
% typically around 40% of the beat.
%
%%%%%

% Set a lower bound of 'ms'
lower_bound = temp_ms;
% Identify the highest local maximum in the second derivative between these two bounds
rel_pks = pks(pks >= lower_bound & pks <= upper_bound);
[~, max_el] = max(curr.derivs.second(rel_pks));
% If this is the first local maximum in this search region ...
if max_el == 1
    % work out whether this has detected the "c" wave
    % - find out if there's an inflection point between "b" and this
    temp_trs = func.waveform.find_pks_trs(curr.derivs.third, 'tr');
    no_infl = sum(temp_trs > temp_b & temp_trs < rel_pks(max_el));
    % - if not then take the next peak
    if no_infl == 0
        % if there is 1 peak in this search region ...
        if length(rel_pks) < max_el+1   % Added in PulseAnalyse5
            % then take the next peak (slightly above the upper bound
            orig_el = find(pks >= lower_bound & pks <= upper_bound);
            try
                rel_pks = pks(orig_el:orig_el+1);
            catch 
                temp_e = [];
                return 
            end
        end
        rel_pk = rel_pks(max_el+1);
    else
        rel_pk = rel_pks(max_el);
    end
else
    rel_pk = rel_pks(max_el);
end
temp_e = rel_pk;

end