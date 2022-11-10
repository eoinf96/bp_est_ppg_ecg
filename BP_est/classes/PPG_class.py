import seaborn as sns
from funcs import HRV, waveform_func
from classes import pw_inds_cls, Waveform
import numpy as np
import matplotlib.pyplot as plt

############################# Class for PPG
class PPG_cls(Waveform):

    def __init__(self,  volunteer_id,study_id='a', device_name='Stowood', folder_loc = None):
        super().__init__(volunteer_id, study_id=study_id,load_mat=0, folder_loc=folder_loc)
        self.onsets = []
        self.peaks = []
        self.ts_filt = []
        self.sqi_beat = []
        self.fid_pts = []
        self.pw_inds = None
        self.norm_fid_pts = [] #Norm features
        self.norm_pw_inds = None
        self.type = 'PPG'
        self.units = 'a.d.u'
        self.device = device_name
        self.study_id = study_id

    def load_mat(self):
        PPG = super().load_mat()
        self._load(PPG=PPG)

    def load_pickle(self, file_name = None):
        PPG = super().load_pickle(file_name=file_name)
        self._load(PPG=PPG)

    def _load(self, PPG):
        self.ts = PPG['ts']
        self.t = PPG['t']
        self.sqi = PPG['sqi']
        self.fs = PPG['fs']
        self.onsets = PPG['onsets']
        self.peaks = PPG['peaks']
        self.ts_filt = PPG['ts_filt']
        self.sqi_beat = PPG['sqi_beat']
        self.fid_pts = PPG['fid_pts']
        self.pw_inds = pw_inds_cls(pw_inds = PPG['pw_inds'], fid_pts=self.fid_pts, sqi_beat = self.sqi_beat, t_beat=self.fid_pts['s']['t'])
        self.norm_fid_pts = PPG['norm_fid_pts']
        self.norm_pw_inds = pw_inds_cls(pw_inds = PPG['norm_pw_inds'], fid_pts=self.norm_fid_pts, sqi_beat = self.sqi_beat, t_beat=self.fid_pts['s']['t'])
        self.t_beat = self.fid_pts['f1']['t']  # Set time as onset time of pulse


    def save_pickle(self, list = None, file_name = None):

        self.pw_inds = return_pw_inds_to_dict(self.pw_inds)
        self.norm_pw_inds = return_pw_inds_to_dict(self.norm_pw_inds)
        if list is None:
            list = ['ts', 't', 'sqi', 'fs', 'onsets', 'peaks', 'fid_pts', 'ts_filt', 'sqi_beat', 'pw_inds', 'norm_fid_pts','norm_pw_inds' ]
        super().save_pickle(list = list, file_name=file_name)

        self.pw_inds = pw_inds_cls(pw_inds = self.pw_inds, fid_pts=self.fid_pts, sqi_beat = self.sqi_beat, t_beat=self.fid_pts['s']['t'])
        self.norm_pw_inds = pw_inds_cls(pw_inds=self.norm_pw_inds, fid_pts=self.fid_pts, sqi_beat=self.sqi_beat, t_beat=self.fid_pts['s']['t'])

    @classmethod
    def get_fid_pts_name(cls):
        return ['g1', 'g2', 'g3', 'g4', 'gauss_error', 'a', 'b', 'c', 'd', 'e', 'f', 's', 'dia', 'dic', 'p1pk', 'p2pk', 'p1in', 'p2in', 'W', 'f1', 'f2', 'halfpoint', 'tangent']


    def _get_fid_pt_ind(self, fid_name):
        return self.fid_pts[fid_name]['ind']

    def get_fid_pt(self, fid_name):
        fid = Waveform(self.volunteer_id, load_mat=0)

        if fid_name in {'g1', 'g2', 'g3', 'g4', 'gauss_error'}:
            print('You need to use get_gauss ')
        else:
            try:
                fid.ts = self.fid_pts[fid_name]['amp']
                fid.t = self.fid_pts[fid_name]['t']
                fid.sqi = np.ones((len(fid.t), 1))
                fid.type =  fid_name
                fid.units = 'unit'
                fid.fs = 0
            except:
                raise Exception('Fid name not known')

        return fid


    ###### PW ind methods

    def get_pw_at_cuff_times(self, BP, pw_names = None, sqi_threshold = 0.8, window_length = 40, do_normalise = True,
                             align_flag = False, return_delay_flag = False, starting_delay = 0):
        initialise_feature_list(self, BP)
        if return_delay_flag:
            print('Havent editted the code yet to return the delay')

        if do_normalise:
            feature_means = self.norm_pw_inds.get_pw_at_cuff_times_no_align(BP, pw_names = pw_names, sqi_threshold = sqi_threshold, window_length = window_length)
            pw_names = self.norm_pw_inds.pw_inds_names
        else:
            feature_means = self.pw_inds.get_pw_at_cuff_times_no_align(BP, pw_names = pw_names, sqi_threshold = sqi_threshold, window_length = window_length)
            pw_names = self.pw_inds.pw_inds_names

        #This is the same order as that processed in the function
        for f_idx, feature_name in enumerate(pw_names):
            self.feature_at_cuff_times[feature_name] = feature_means[:,f_idx]
            self.feature_BP_correlation[feature_name] = waveform_func.get_correlation_of_signals(self.feature_at_cuff_times[feature_name], BP.ts)


    @classmethod
    def get_pw_inds_name(cls):
        return pw_inds_cls.get_pw_inds_name()
    @classmethod
    def get_gauss_pw_inds_name(cls):
        return pw_inds_cls.get_gauss_pw_inds_name()

    @classmethod
    def get_Monte_Moreno_feature_names(cls):
        feature_name_list = waveform_func.get_Monte_Moreno_features(None, return_feature_list=True)
        return feature_name_list
    @classmethod
    def get_HRV_features(cls):
        return ['PPG_SDRR', 'PPG_pRR50', 'PPG_RMSSD', 'PPG_nLF_power', 'PPG_nHF_power',
                'PPG_T_power', 'PPG_LF_div_HF']
    @classmethod
    def get_windowed_features_list(cls, do_MM = 1, do_HRV = 1):
        feature_name_list = []
        if do_MM : feature_name_list = feature_name_list + cls.get_Monte_Moreno_feature_names()
        if do_HRV: feature_name_list = feature_name_list + cls.get_HRV_features()
        return feature_name_list

    def get_windowed_features(self, BP, window_length = 60, good_beat_threshold = 0.5, sqi_threshold = 0.8, feature_names ='All'):
        do_MM = 0
        do_HRV = 0
        if feature_names.upper() == 'ALL':
            do_MM = 1
            do_HRV = 1
        elif feature_names.upper() == 'MONTE_MORENO':
            do_MM = 1
        elif feature_names.upper() == 'HRV':
            do_HRV = 1
        initialise_feature_list(self, BP)

        data = []
        window_start_times = BP.t[BP.sqi > 0] - window_length / 2
        window_end_times = BP.t[BP.sqi > 0] + window_length / 2

        if do_MM:
            features = local_func_get_MM_at_cuff_times(self, window_start_times= window_start_times, window_end_times = window_end_times,good_beat_threshold = good_beat_threshold, sqi_threshold = sqi_threshold)
            data.extend(features.values())

        if do_HRV:
            t_point = self.fid_pts['tangent']['t']
            t_point = t_point[~np.isnan(t_point)]

            RR_t = t_point[1:]
            RR_ts = t_point[1:] - t_point[:-1]
            features = HRV.get_HRV_features(RR_ts_beat=RR_ts, RR_t_beat=RR_t, window_start_times= window_start_times, window_end_times=window_end_times)

            data.extend(features.values())

        data = np.vstack(data)
        #Get inflation times with error in them
        no_error = False
        try:
            inf_error = (data == 'Error').any(axis =0)
        except AttributeError:
            no_error = True
        if not no_error:
            median_values = np.nanmedian(data[:, ~inf_error].astype('float'), axis = 1)
            data[:, inf_error] = np.array([median_values]*sum(inf_error)).T
            data = data.astype('float')
        for f_idx, feat_name in enumerate(self.get_windowed_features_list(do_MM=do_MM, do_HRV=do_HRV)):
            self.feature_at_cuff_times[feat_name] = data[f_idx,:]
            self.feature_BP_correlation[feat_name] = waveform_func.get_correlation_of_signals(data[f_idx,:], BP.ts)




    #Override parent function for segment data for specific features of PPG
    def segment_data(self, t_start, t_end, do_fid_pts = True, do_pw_inds = True, do_peaks = True):
        Out_seg = self.copy()

        # Adjust the t_start and t_end
        t_start = max(t_start, self.t[0])
        t_end = min(t_end, self.t[-1])

        if t_start > t_end:
            Out_seg.t = []
            Out_seg.ts = []
            Out_seg.sqi = []

            return Out_seg

        # This function segments the waveform time series at the input times
        index_ts = np.where(np.logical_and(self.t >= t_start, self.t <= t_end))
        if len(index_ts) == 0:
            raise Exception('t_start and t_end are not in the signal time vector')

        Out_seg.t = Out_seg.t[index_ts]
        Out_seg.t = Out_seg.t - Out_seg.t[0]
        Out_seg.ts = Out_seg.ts[index_ts]
        Out_seg.sqi = Out_seg.sqi[index_ts]

        if do_peaks and hasattr(self, 'peaks') and hasattr(self, 'onsets'):
            index_peak = np.where(np.logical_and(self.peaks >= index_ts[0][0], self.peaks <= index_ts[0][-1]))
            index_peak = index_peak[0]
            index_onsets = np.where(np.logical_and(self.onsets >= index_ts[0][0], self.onsets <= index_ts[0][-1]))
            index_onsets = index_onsets[0]

            if self.peaks[index_peak[0]] < self.onsets[index_onsets[0]]:
                index_peak = np.delete(index_peak, 0)

            if self.peaks[index_peak[-1]] > self.onsets[index_onsets[-1]]:
                index_peak = np.delete(index_peak, -1)
            Out_seg.peaks = Out_seg.peaks[index_peak]
            Out_seg.peaks -= index_ts[0][0]
            Out_seg.onsets = Out_seg.onsets[index_onsets]
            Out_seg.onsets -= index_ts[0][0]
            Out_seg.sqi_beat =Out_seg.sqi_beat[index_peak]


        if do_fid_pts:
            if hasattr(self, 'fid_pts'):
                fid_names = list(self.fid_pts.keys())
                for gauss_name in ['g1', 'g2', 'g3', 'g4', 'gauss_error']:
                    if gauss_name in fid_names: fid_names.remove(gauss_name)

                for fid_name in fid_names:
                    Out_seg.fid_pts[fid_name]['ind'] = Out_seg.fid_pts[fid_name]['ind'][index_peak]
                    Out_seg.fid_pts[fid_name]['amp'] = Out_seg.fid_pts[fid_name]['amp'][index_peak]
                    Out_seg.fid_pts[fid_name]['t'] = Out_seg.fid_pts[fid_name]['t'][index_peak] - t_start


        if do_pw_inds:
            if hasattr(self, 'pw_inds'):
                pw_names = list(self.pw_inds.keys())

                for pw_name in pw_names:
                    Out_seg.pw_inds[pw_name] = Out_seg.pw_inds[pw_name][index_peak]

            if 'g1' in self.fid_pts:
                Out_seg.fid_pts['gauss_error'] = Out_seg.fid_pts['gauss_error'][index_peak]

                for gauss_num in ['g1', 'g2', 'g3']:
                    for var_name in ['amp', 'mu', 'sigma']:
                        exec('Out_seg.fid_pts["'+gauss_num+'"]["'+var_name+'"]    = Out_seg.fid_pts["'+gauss_num+'"]["'+var_name+'"][index_peak]')

        return Out_seg


    def show(self, mins_flag = True, fid_point='f1'):

        fid_loc = self.fid_pts[fid_point]['ind']
        fid_loc = fid_loc[~np.isnan(fid_loc)]
        factor = 60 if mins_flag else 1

        sns.set_style("darkgrid")
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax1.plot(self.t/factor, self.ts, zorder = 1)
        ax1.scatter(np.take(self.t, fid_loc.astype('int64'))/factor, np.take(self.ts, fid_loc.astype('int64')), color='red', zorder =2, label=fid_point)
        plt.ylabel('PPG (a.d.u)', fontsize=16)

        ax2.plot(self.t/factor, self.sqi,zorder = 1)

        xlabel_str = 'Time (mins)' if mins_flag else 'Time (s)'

        plt.xlabel(xlabel_str, fontsize=12)
        plt.ylabel('SQI', fontsize=12)
        plt.legend()


        plt.show()


    def show_pw_indices(self):
        for name in self.pw_inds.keys():
            print(name)



def initialise_feature_list(PPG, BP):
    if not hasattr(PPG, 'feature_at_cuff_times'):
        PPG.feature_at_cuff_times = {
            BP.type: BP.ts
        }
        PPG.feature_BP_correlation = {}


def local_func_get_MM_at_cuff_times(PPG, window_start_times, window_end_times,good_beat_threshold = 0.5, sqi_threshold = 0.8):

    Monte_moreno_features_at_cuff_times = {feature: [] for feature in PPG.get_Monte_Moreno_feature_names()}

    for wind_idx, wind_start in enumerate(window_start_times):
        if wind_idx == 26:
            a=1
        wind_end = min(PPG.t[-1], window_end_times[wind_idx])
        wind_start = max(0, wind_start)

        if wind_end - wind_start < 3:
            for feature in PPG.get_Monte_Moreno_feature_names():
                Monte_moreno_features_at_cuff_times[feature].append('Error')
            continue
        signal_seg = PPG.segment_data(t_start=wind_start, t_end=wind_end, do_pw_inds=False, do_fid_pts=False)
        if len(signal_seg.ts) == 0:
            for feature in PPG.get_Monte_Moreno_feature_names():
                Monte_moreno_features_at_cuff_times[feature].append('Error')
            continue
        # Check if more than good_beat_threshold proportion are of good quality
        if len(signal_seg.get_good_values()) / len(signal_seg.ts) < good_beat_threshold:
            for feature in PPG.get_Monte_Moreno_feature_names():
                # Monte_moreno_features_at_cuff_times[feature].append('Error')
                Monte_moreno_features_at_cuff_times[feature].append(np.nan)
            continue

        features_window = waveform_func.get_Monte_Moreno_features(signal_seg)
        for feature in PPG.get_Monte_Moreno_feature_names():
            Monte_moreno_features_at_cuff_times[feature].append(features_window[feature])

    return Monte_moreno_features_at_cuff_times

#Function needed in order to properly save pw_inds
def return_pw_inds_to_dict(pw_inds):
    out_dict = {}
    for pw_ind_name in pw_inds.pw_inds_names:
        exec('out_dict[pw_ind_name] = pw_inds.'+pw_ind_name)
    return out_dict