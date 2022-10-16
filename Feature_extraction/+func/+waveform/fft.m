function [pw, pw_freq] = fft(sig_in, fs, opts)
narginchk(2, inf)
%% Update opts
if nargin < 3
    opts  = struct();
end
default_opts.window_data = 0;
default_opts.detrend =1;
opts = func.aux_functions.update_with_default_opts(opts, default_opts);
%% Preprocess
sig_in = sig_in(:);	
N_sig = length(sig_in);
if opts.detrend 
    sig_in =detrend(sig_in, 0, 'omitnan');
end
if opts.window_data
    window_data = hamming( N_sig );
    sig_in = sig_in .* window_data;
end
%% Compute FFT
nfft = 2^nextpow2( N_sig );  
y = fft(sig_in, nfft);

% frequency vector
pw_freq = (0:nfft/2)' * fs/nfft; %f(1:N);

% Take the magnitude of fft and scale it
pw_mag = abs( y(1:(nfft/2+1)) / nfft );
sq_mag = pw_mag .^2;
pw = sq_mag / norm(sq_mag);
end

