from classes import Waveform, PPG_cls, ECG_cls
import seaborn as sns
from funcs import waveform_func
import numpy as np
import matplotlib.pyplot as plt
#Allow for loading
#And allow for complete recomputing via ECG or PPG if neccesary
class PAT_cls(Waveform):
    def __init__(self,  volunteer_id, device_name='Stowood', study_id = 'a'):
        super().__init__(volunteer_id, load_mat=0 ,  study_id = study_id)
        self.t_beat = []
        self.ts_beat = []
        self.sqi_beat = []
        self.type = 'PAT'
        self.units = 's'
        self.distal_point = None
        self.device = device_name
        self.study_id = study_id

    def load_mat(self):
        PAT = super().load_mat()
        self._load(PAT=PAT)

    def load_pickle(self, file_name = None):
        PAT = super().load_pickle(file_name=file_name)
        self._load(PAT=PAT)

    def _load(self, PAT):
        self.ts = PAT['ts']
        self.t = PAT['t']
        self.sqi = PAT['sqi']
        self.fs = PAT['fs']
        self.sqi_beat = PAT['sqi_beat']
        self.ts_beat = PAT['ts_beat']
        self.t_beat = PAT['t_beat']

    def save_pickle(self,list = None, file_name = None):
        if list is None:
            list = ['ts', 't', 'sqi', 'fs', 'sqi_beat', 'ts_beat', 't_beat']
        super().save_pickle(list=list, file_name=file_name)



    def get_PAT_beat(self, PPG=None, ECG=None, distal_name = 'tangent'):
        if PPG is None:
            #Then we need to load the PPG
            PPG = PPG_cls(volunteer_id = self.volunteer_id)
            PPG.load_mat()
        if ECG is None:
            ECG = ECG_cls(volunteer_id = self.volunteer_id)
            ECG.load_mat()

        #Loop through all ECG peaks and find the PPG distal point
        distal_points = PPG.get_fid_pt(distal_name)
        num_beats = len(ECG.peaks) -1

        #initalise variables
        self.t_beat = ECG.index_t(ECG.peaks[:-1]).reshape((num_beats,))
        self.ts_beat = np.empty((num_beats, 1))
        self.ts_beat[:] = np.nan
        self.sqi_beat = np.ones((num_beats, 1))


        for ecg_idx in range(num_beats):

            if ECG.sqi_beat[ecg_idx]==0:
                continue
            ecg_start_time  = ECG.index_t(ECG.peaks[ecg_idx])
            ecg_end_time = ECG.index_t(ECG.peaks[ecg_idx+1])

            ppg_beat_index = np.argwhere(np.logical_and(distal_points.t >ecg_start_time,distal_points.t < ecg_end_time))
            if ppg_beat_index.size == 0:
                self.sqi_beat[ecg_idx] = 0
            else:
                ppg_beat_index = ppg_beat_index[0]
                self.sqi_beat[ecg_idx] = min(ECG.sqi_beat[ecg_idx], PPG.sqi_beat[ppg_beat_index])
                self.ts_beat[ecg_idx] = distal_points.t[ppg_beat_index] - ecg_start_time

        #remove any beats that still have nan as their ts
        self.t_beat = self.t_beat[~np.isnan(self.ts_beat).reshape(num_beats,)]
        self.sqi_beat = self.sqi_beat[~np.isnan(self.ts_beat)]
        self.ts_beat = self.ts_beat[~np.isnan(self.ts_beat)]

    def process_PAT_beats(self, sqi_threshold = 0.8, debug = False):
        #Run process from the func file
        waveform_func.process_beat_data_MOLLIE(self, sqi_threshold = sqi_threshold, debug=debug)

    #Overwrite of parent waveform function so that we dont have to state to inverse PAT
    def get_signal_at_cuff_times(self, BP, align_flag = True, invert_flag = True,  max_lag = 5*60,
                                 window_size = 30, sqi_threshold = 0.8,
                                 good_beat_threshold = 0.5, debug = False, debug_deep = False,
                                 return_correlation_flag = False, return_delay_flag = False):

        invert_flag = True
        output = super().get_signal_at_cuff_times(BP, align_flag, invert_flag,  max_lag,
                                 window_size, sqi_threshold,
                                 good_beat_threshold, debug, debug_deep ,
                                 return_correlation_flag, return_delay_flag)
        return output



    def show(self, mins_flag=True):
        factor = 60 if mins_flag else 1

        sns.set_style("darkgrid")
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ind_start = 4 * 60 * 60
        ax1.plot(self.t[ind_start:] / factor, self.ts[ind_start:], zorder=1)
        # ax1.plot(self.t / factor, self.ts, zorder=1)
        plt.ylabel('PAT (s)', fontsize=16)

        ax2.plot(self.t[ind_start:] / factor, self.sqi[ind_start:], zorder=1)

        xlabel_str = 'Time (mins)' if mins_flag else 'Time (s)'

        plt.xlabel(xlabel_str, fontsize=12)
        plt.ylabel('SQI', fontsize=12)

    def show_beat(self, mins_flag=True):
        factor = 60 if mins_flag else 1

        sns.set_style("darkgrid")
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax1.plot(self.t_beat / factor, self.ts_beat, zorder=1)
        # ax1.plot(self.t / factor, self.ts, zorder=1)
        plt.ylabel('PAT (s)', fontsize=16)

        ax2.plot(self.t_beat / factor, self.sqi_beat, zorder=1)

        xlabel_str = 'Time (mins)' if mins_flag else 'Time (s)'

        plt.xlabel(xlabel_str, fontsize=12)
        plt.ylabel('SQI', fontsize=12)




    def print_hello(self):
        print('Hello')
