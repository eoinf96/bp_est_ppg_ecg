function temp_p1 = identify_p1(curr, temp_b, buffer_p1)

% find p1

% find local minima in the first derivative
fd_trs = func.waveform.find_pks_trs(curr.derivs.first, 'tr');


% Find the first local minimum in the first derivative after 0.1 s
current_buffer = buffer_p1(1);
temp = find(fd_trs > current_buffer,1);
% Find the second local minimum (or failing that, the first) in the first derivative after 'b'
temp2 = find(fd_trs > temp_b, 2);
if length(temp2) > 1
    temp2 = temp2(2);
end
% Take whichever comes first:
if temp2 < temp
    temp = temp2;
end
temp_p1 = fd_trs(temp);

% If this value for p1 is after the buffer of 0.18 s ...
if temp_p1 > buffer_p1(2)
    
    % Then find the first local minimum in the fourth derivative after 0.1 s
    fd_trs = func.waveform.find_pks_trs(curr.derivs.fourth, 'tr');
    temp = find(fd_trs > current_buffer,1);
    temp_p1 = fd_trs(temp); clear temp
end

% If this value for p1 is after the buffer of 0.18 s ...
if temp_p1 > buffer_p1(2)
    % Then find the last local minimum in the first derivative before 0.18 s
    temp_p1 = fd_trs(find(fd_trs <= current_buffer,1,'last'));
   
    % If this doesn't find temp_p1, then extend the buffer
    if isempty(temp_p1)
        temp_p1 = fd_trs(find(fd_trs > current_buffer,1,'first'));
    end
end



end