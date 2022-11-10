from classes import Waveform
from funcs import waveform_func
import numpy as np

class pw_inds_cls:
    def __init__(self, pw_inds, fid_pts, sqi_beat, t_beat):
        # self.pw_inds = pw_inds
        self.fid_pts = fid_pts
        self.sqi_beat = sqi_beat
        self.num_pw_inds = len(pw_inds)
        self.t_beat =t_beat

        self.pw_inds_names = []
        for key in pw_inds:
            exec('self.'+key+'= pw_inds[key]')
            self.pw_inds_names += [key]
        self.pw_inds_names = [x for x in self.get_pw_inds_name() if x in self.pw_inds_names]
        # self.pw_inds_names_missed = [x for x in self.get_pw_inds_name() if x not in self.pw_inds_names]


    def __getitem__(self, item):
        return getattr(self, item )

    def update_pw_inds_names(self, added_features):
        self.pw_inds_names +=added_features
        self.pw_inds_names = list(set(self.pw_inds_names))
        self.pw_inds_names = [x for x in self.get_pw_inds_name() if x in self.pw_inds_names]
        # self.pw_inds_names_missed = [x for x in self.get_pw_inds_name() if x not in self.pw_inds_names]
        self.num_pw_inds = len(self.pw_inds_names)


    def get_pw_ind(self, pw_name, process_flag = False, units = 'NA', plot_flag = False):

        pw_ind = Waveform(volunteer_id=0, load_mat=0)

        pw_ind.type = pw_name
        pw_ind.units = units
        pw_ind.fs = 0
        pw_ind.ts_beat = self[pw_name]
        pw_ind.t_beat = self.t_beat  # Set time as onset time of pulse

        if pw_name in self.get_gauss_pw_inds_name():
            #Then we set SQI by GE threshold
            pw_ind.sqi_beat = self.sqi_beat * get_gauss_SQI(self.fid_pts['gauss_error'], error_threshold = 0.023)
        else:
            pw_ind.sqi_beat = self.sqi_beat

        #Remove any nans
        isnan = np.isnan(pw_ind.ts_beat)
        med = np.median(pw_ind.ts_beat[~isnan])
        pw_ind.sqi_beat[isnan] = 0
        pw_ind.ts_beat[isnan] = med

        if process_flag:
            #from Func
            pw_ind = waveform_func.process_beat_data_MOLLIE(pw_ind, debug = plot_flag)

        return pw_ind



    def get_pw_at_cuff_times(self, BP, pw_names = None, sqi_threshold = 0.8,
                             align_flag = False, return_delay_flag = False, starting_delay = 0):
        #Define a new property that contains the pw inds at cuff value
        feature_at_cuff_times = {
            BP.type: BP.ts
        }
        feature_BP_correlation = {}

        if return_delay_flag: feature_delay = {}

        if pw_names is None:
            pw_names = self.get_pw_inds_name()
        if pw_names is 'gauss':
            #Then do Gauss PW only
            pw_names = self.get_gauss_pw_inds_name()
        self.pw_inds_names = pw_names
        if align_flag:
            if not hasattr(BP, 'align'):
                BP.get_BP_align()
            BP_align = BP.align

        for p_idx, pw_name in enumerate(pw_names):
            print(pw_name)
            pw = self.get_pw_ind(pw_name, process_flag = True)

            if len(pw.ts)==0:
                #then the processed signal has no ts so set correlation to nan
                sig_cuff = []
                correlation = np.nan
                delay=0
            else:
                delay = starting_delay
                if align_flag:
                    #Look into the invert flag
                    delay = waveform_func.find_delay_signals(Signal_one=BP_align, Signal_two=pw, invert_two_flag=0)
                    delay = delay/pw.fs
                pw.t = pw.t + delay

                #Get pw values at times of cuff inflation
                sig_cuff, correlation = pw.get_signal_at_cuff_times(BP = BP, align_flag=False, sqi_threshold=sqi_threshold, return_correlation_flag = True)
                sig_cuff = np.vstack(sig_cuff)
                sig_cuff = sig_cuff[:,1]

            feature_at_cuff_times[pw_name] = sig_cuff
            feature_BP_correlation[pw_name] = correlation
            if return_delay_flag: feature_delay[pw_name] = delay

        if not return_delay_flag:
            return feature_at_cuff_times, feature_BP_correlation
        else:
            return feature_at_cuff_times, feature_BP_correlation, feature_delay

    def get_pw_at_cuff_times_no_align(self, BP, pw_names=None, sqi_threshold=0.8, window_length = 40):

        t_start_times = BP.t - window_length/2
        t_end_times = BP.t + window_length / 2

        #Change this to be the number currently loaded??
        if pw_names is None:
            #Then do all
            pw_names = self.get_pw_inds_name()
        if pw_names is 'gauss':
            #Then do Gauss PW only
            pw_names = self.get_gauss_pw_inds_name()
        self.pw_inds_names = pw_names
        feature_stack = self[pw_names[0]]
        for pw_name in pw_names[1:]:
            pw = self[pw_name]
            if pw_name in self.get_gauss_pw_inds_name():
                sqi_beat = self.sqi_beat * get_gauss_SQI(self.fid_pts['gauss_error'], error_threshold=0.023)
            else:
                sqi_beat = self.sqi_beat
            pw[sqi_beat < sqi_threshold] = np.nan
            feature_stack = np.vstack((feature_stack, self[pw_name]))

        median_feature_vals = np.nanmedian(feature_stack, axis = 1)
        mean_store = []
        for idx, t_start in enumerate(t_start_times):
            t_end = t_end_times[idx]
            if t_end < self.t_beat[~np.isnan(self.t_beat)][0] or t_start > self.t_beat[~np.isnan(self.t_beat)][-1]:
                mean_store.append(median_feature_vals)
                continue

            ind_in_window = np.where(np.logical_and(self.t_beat > t_start, self.t_beat < t_end))

            mean_store.append(np.squeeze(np.nanmean(feature_stack[:, ind_in_window], axis = 2)))


        feature_means =np.array(mean_store)

        return feature_means


    #This is just a guess and may need to chabge depending on if the pw inds have dynamically changed
    @classmethod
    def get_pw_inds_name(cls):
        timing_features = ['delta_t', 'CT', 't_sys', 't_dia', 't_ratio']
        amplitude_features = ['dic_amp', 'RI', 'sVRI']
        area_width_features = ['A1', 'A2', 'IPA', 'width25', 'width50']
        first_deriv_features = ['sp_mean', 'sp_var', 'dp_mean', 'dp_var']
        # second_deriv_features = ['b_div_a', 'c_div_a', 'd_div_a', 'e_div_a', 'a_div_amp', 'b_div_amp', 'c_div_amp',
        #                          'd_div_amp', 'e_div_amp', 'AGI', 'slope_b_c', 'slope_b_d', 'PPG_AI']
        second_deriv_features = ['b_div_a', 'c_div_a', 'd_div_a', 'e_div_a','AGI', 'slope_b_c', 'slope_b_d', 'PPG_AI']
        mix_deriv_features = ['PI', 'STT']
        freq_features = ['NHA', 'skewness', 'kurtosis', 'IHAR']
        # liang_features = ['liang_1', 'liang_2', 'liang_3', 'liang_4', 'liang_5','liang_6', 'liang_7', 'liang_8', 'liang_9', 'liang_10']
        pw_inds = []
        pw_inds += timing_features
        pw_inds += amplitude_features
        pw_inds += area_width_features
        pw_inds += first_deriv_features
        pw_inds += second_deriv_features
        pw_inds += mix_deriv_features
        pw_inds += freq_features
        # pw_inds += liang_features
        pw_inds += cls.get_gauss_pw_inds_name()

        return pw_inds

    @classmethod
    def get_gauss_pw_inds_name(cls):
        return ['gauss_AI', 'gauss_RI', 'gauss_RTT_Rubins', 'gauss_AIx_Rubins', 'gauss_RI_Rubins', 'gauss_LVET', 'gauss_sys_dias', 'gauss_amp4_amp1', 'gauss_sigma4_amp1', 'g1_amp', 'g1_sigma', 'g1_mu', 'g2_amp', 'g2_sigma', 'g2_mu', 'g3_amp', 'g3_sigma', 'g3_mu', 'g4_amp', 'g4_sigma', 'g4_mu']



def get_gauss_SQI(GE, error_threshold = 0.023):

    sqi_beat = np.ones((len(GE),))
    sqi_beat[np.where(GE > error_threshold)] = 0
    return sqi_beat