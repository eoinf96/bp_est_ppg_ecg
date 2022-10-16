function derivs = get_derivs(ts, fs)
s_g_filter_len = 7;
dt = 1/fs;
derivs.first  = func.waveform.savitzky_golay_deriv(ts, 1, s_g_filter_len)./dt;
derivs.second = func.waveform.savitzky_golay_deriv(ts, 2, s_g_filter_len)./dt./dt;
derivs.third  = func.waveform.savitzky_golay_deriv(ts, 3, s_g_filter_len)./dt./dt./dt;
derivs.fourth = func.waveform.savitzky_golay_deriv(ts, 4, s_g_filter_len)./dt./dt./dt./dt;
end

