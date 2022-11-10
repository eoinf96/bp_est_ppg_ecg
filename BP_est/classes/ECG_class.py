import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from funcs import filtering_funcs, HRV
import matlab
import matlab.engine
from classes import Waveform


############################# Class for ECG
class ECG_cls(Waveform):
    def __init__(self,  volunteer_id, device_name='Stowood', study_id = 'a'):

        super().__init__(volunteer_id, load_mat=0, study_id =study_id)
        self.onsets = []
        self.peaks = []
        self.fid_pts = []
        self.ts_filt = []
        self.sqi_beat = []
        self.pw_inds = []
        self.type = 'ECG'
        self.units = 'mV'
        self.device = device_name
        self.study_id = study_id

    def load_mat(self):
        ECG = super().load_mat()
        self._load(ECG=ECG)

    def load_pickle(self, file_name = None):
        ECG = super().load_pickle(file_name=file_name)
        self._load(ECG=ECG)

    def _load(self, ECG):
        #z normalise the ECG -- inorder to normalise the y axis from the 3 devices
        ecg = (ECG['ts'] - np.mean(ECG['ts']))/(np.std(ECG['ts']))
        # ecg = ECG['ts']

        # self.ts = ECG['ts']
        self.ts = ecg
        self.t = ECG['t']
        self.sqi = ECG['sqi']
        self.fs = ECG['fs']
        self.onsets = ECG['onsets']
        self.peaks = ECG['peaks']
        self.sqi_beat = ECG['sqi_beat']

    def save_pickle(self, list = None, file_name = None):
        if list is None:
            list = ['ts', 't', 'sqi', 'fs', 'onsets', 'peaks', 'sqi_beat', ]
        super().save_pickle(list=list, file_name=file_name)

    def get_features_at_cuff_times(self, BP, window_size=40, sqi_value = 0.8, debug = False, feature_names = 'All', eng = None):
        do_sim = 0
        do_HRV = 0

        if isinstance(feature_names, str):
            if feature_names.upper() == 'ALL':
                do_sim = 1
                do_HRV = 1
                feature_names = None
            elif feature_names.upper() == 'SIM':
                do_sim = 1
                feature_names = None
            elif feature_names.upper() == 'HRV':
                do_HRV = 1
                feature_names = None

        else:
            if len(list(set(feature_names) & set(self.get_features_list(do_HRV = 0)))) > 0:
                do_sim = 1
            elif len(list(set(feature_names) & set(self.get_features_list(do_sim = 0)))) > 0:
                do_HRV = 1


        self.feature_at_cuff_times = {
            BP.type: BP.ts
        }
        # self.feature_BP_correlation = {}

        data = []

        #Do windows around the BP values that are not sqi = 0
        window_start_times = BP.t[BP.sqi > 0] - window_size / 2
        window_end_times = BP.t[BP.sqi > 0] + window_size / 2

        if do_sim:
            features = get_simjanoska_features(self, window_start_times= window_start_times, window_end_times=window_end_times, feature_list = feature_names, eng = eng )
            data.extend(features.values())

        if do_HRV:
            RR_t = self.t[self.peaks[1:]]
            RR_ts = self.t[self.peaks[1:]] - self.t[self.peaks[:-1]]
            features = HRV.get_HRV_features(RR_ts_beat=RR_ts, RR_t_beat=RR_t, window_start_times= window_start_times, window_end_times=window_end_times)#, feature_list = feature_names )
            data.extend(features.values())

        data = np.vstack(data)
        if feature_names is None: feature_names = self.get_features_list(do_sim=do_sim, do_HRV=do_HRV)
        if debug:
            fig, ax = plt.subplots(len(feature_names) +1,1)
            ax[0].plot(BP.t, BP.ts)
            ax[0].set_ylabel('BP mmHg')
            ax[0].set_xlabel('Time (s)')

            for idx in range(len(feature_names)):
                axes = ax[idx+1]
                axes.plot(BP.t, data[:,idx])
                axes.set_ylabel(feature_names[idx])
                axes.set_xlabel('Time (s)')

        for f_idx, feat_name in enumerate(feature_names):
            self.feature_at_cuff_times[feat_name] = data[f_idx,:]
            # self.feature_BP_correlation[feat_name] = get_correlation_of_signals(data[f_idx,:], BP.ts)


    def remove_power_line_interference(self):
        fs_power_line = 50
        self.ts= filtering_funcs.notch(self.ts, fs_power_line, self.fs, q=35)

    def remove_baseline(self):
        #For now using suggestions from rubbish paper on BP estimation with ECG -- Simjanoska

        fs_cut_off = 0.3
        order = 4
        self.ts = filtering_funcs.highpass(self.ts, fs = self.fs,cutoff_high =fs_cut_off,  order=order)

    @classmethod
    def get_features_list(cls, do_sim = 1, do_HRV = 1):
        feature_list = []
        if do_sim: feature_list +=['Complexity', 'Mobility', 'Shannon', 'Fractal', 'Holder', 'Sec_cummulant']
        if do_HRV: feature_list +=['SDRR', 'pRR50', 'RMSSD', 'nLF_power',  'nHF_power', 'T_power' , 'LF_div_HF']
        return feature_list


        # return complexity, mobility, Shannon, fractal

    def show(self, mins_flag = True):
        factor = 60 if mins_flag else 1

        sns.set_style("darkgrid")
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax1.plot(self.t / factor, self.ts, zorder=1)
        ax1.scatter(self.index_t(self.peaks) / factor, self.index_ts(self.peaks), color='red', zorder=2)
        plt.ylabel('ECG (mV)', fontsize=16)

        ax2.plot(self.t / factor, self.sqi, zorder=1)

        xlabel_str = 'Time (mins)' if mins_flag else 'Time (s)'

        plt.xlabel(xlabel_str, fontsize=12)
        plt.ylabel('SQI', fontsize=12)

        plt.show()



def hFD(a, k_max = 282): #Higuchi FD

    L = []
    x = []
    N = len(a)

    for k in range(1,k_max):
        Lk = 0
        for m in range(0,k):
            #we pregenerate all idxs
            idxs = np.arange(1,int(np.floor((N-m)/k)),dtype=np.int32)
            Lmk = np.sum(np.abs(a[m+idxs*k] - a[m+k*(idxs-1)]))
            Lmk = (Lmk*(N - 1)/(((N - m)/ k)* k)) / k
            Lk += Lmk

        L.append(np.log(Lk/(m+1)))
        x.append([np.log(1.0/ k), 1])

    (p, r1, r2, s)=np.linalg.lstsq(x, L, rcond=None)
    return p[0]

def get_simjanoska_features(ECG, window_start_times, window_end_times, feature_list = None, good_beat_threshold=0.5, eng = None):


    if feature_list is None: feature_list = ['Complexity', 'Mobility', 'Shannon', 'Fractal', 'Holder', 'Sec_cummulant']
    #Get lowercase of feature_list
    feature_list = [x.lower() for x in feature_list]

    result = {feature: [] for feature in feature_list}

    for wind_idx, wind_start in enumerate(window_start_times):
        wind_end = min(ECG.t[-1], window_end_times[wind_idx])
        wind_start = max(0, wind_start)

        ECG_seg = ECG.segment_data(t_start=wind_start, t_end=wind_end)

        if len(ECG_seg.ts) == 0:
            #Then the ECG does not overlap with this cuff inflation
            for feature in feature_list:
                # result[feature].append('Error')
                result[feature].append(np.nan)
            continue
        # If the data isnt of good enough quality then we will have to linearly interpolate at the end
        if len(ECG_seg.get_good_values()) / len(ECG_seg.ts) < good_beat_threshold:
            for feature in feature_list:
                result[feature].append(np.nan)
            continue

        #make the signal start and end at the last part of good quality. Any bad quality in between nothing can be done sadly
        good_idx = [x for x,v in enumerate(ECG_seg.sqi) if v>0.8]
        start_idx = good_idx[0]
        end_idx = good_idx[-1]
        ECG_seg.ts=ECG_seg.ts[start_idx:end_idx]


        if 'complexity' in feature_list or 'mobility' in feature_list:
            dts = np.diff(np.append(0, ECG_seg.ts))
            ddts = np.diff(np.append(0, dts))

            mts2 = np.mean(ECG_seg.ts **2)
            mdts2 = np.mean(dts ** 2)
            mddts2 = np.mean(ddts ** 2)

            mob = mdts2 /mts2

            result['complexity'].append(np.sqrt((mddts2/ mdts2) - mob))
            result['mobility'].append(np.sqrt(mob))

        if 'shannon' in feature_list:
            probs, bin_edges = np.histogram(ECG_seg.ts, density=True)
            # Add correction term for bin width
            width = np.median(np.diff(bin_edges)) # Should all be the same but take median incase
            result['shannon'].append(-np.sum((probs * np.log2(probs)) + width))

        if 'fractal' in feature_list: result['fractal'].append(hFD(ECG_seg.ts, k_max = 282))

        if 'holder' in feature_list:
            ts = matlab.double(ECG_seg.ts.tolist())
            dh, h, cp = eng.dwtleader(ts, nargout=3)
            result['holder'].append(h[0][np.argmax(dh)])
            result['sec_cummulant'].append(cp[0][1])


    return result