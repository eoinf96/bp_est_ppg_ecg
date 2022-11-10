from loading_functions import *
import seaborn as sns
import copy
from set_root import set_root
import scipy.fftpack
import os
import _pickle as cPickle
import pickle


class Waveform:
    root = set_root()

    def __init__(self, volunteer_id, study_id = 'a', load_mat = 1, units = [], folder_loc = None):
        if not type(volunteer_id) == str:
            volunteer_id = '00' + str(volunteer_id)
            volunteer_id = volunteer_id[-3:]

        self.volunteer_id = volunteer_id
        self.study_id = study_id
        if folder_loc is None:
            self.file_loc = Waveform.root + '/Studies/MOLLIE/processed_signals/' + self.volunteer_id + self.study_id + '.mat'
        else:
            self.file_loc = Waveform.root + folder_loc + self.volunteer_id + self.study_id + '.mat'
        self.record_name = self.volunteer_id + self.study_id
        self.ts = []
        self.t = []
        self.sqi = []
        self.fs = []
        self.units = units
        self.type = ''
        self.device = ''


        if load_mat:
            self._mat_contents = sio.loadmat(self.file_loc, squeeze_me=True)

    def load_mat(self):
        self.record_name = self.volunteer_id + self.study_id
        print("Loading "+ self.type + ": " + self.file_loc)
        self._mat_contents = sio.loadmat(self.file_loc, squeeze_me=True)
        #Change get -- to include any signal
        sig = 0
        if self.type in ['SBP', 'MAP', 'DBP']:
            sig = get_BP(self._mat_contents, self.type)
        else:
            sig = eval('get_'+self.type+'(self._mat_contents, self.device)')
        return sig

    def load_pickle(self, file_name = None):
        if not self.type[1:] == 'BP':
            if self.device.upper() =='STOWOOD':
                ending = self.type+'_s.dat'
            elif self.device.upper() == 'PHILIPS':
                ending = self.type+'_p.dat'
            elif self.device.upper() == 'CARDIOSCREEN':
                ending = self.type+'_c.dat'
            else:
                ValueError('Unknown device name')

        else:
            ending = self.type+'.dat'


        if file_name is None:
            self.file_loc = Waveform.root + '/Studies/MOLLIE/processed_signals/' + self.volunteer_id + self.study_id + '_pickle/' + ending
        else:
            self.file_loc = Waveform.root + '/Studies/MOLLIE/processed_signals/' + file_name + '_pickle/' + ending
        print("Loading " + self.type + ": " + self.file_loc)
        # Method to load the BP class as a pickle file rather than mat file
        if os.path.getsize(self.file_loc) > 0:
            with open(self.file_loc, 'rb') as inpf:
                sig = cPickle.load(inpf)
        else:
            raise ValueError('No file saved under the filename')
        return sig

    def save_pickle(self, list, file_name = None):
        # Method to save the PPG class as a pickle file -- have not saved as class instance as will not be needing to save the methods etc
        # Be careful though!! Can override data
        if self.device.upper() == 'STOWOOD':
            ending = self.type + '_s.dat'
        elif self.device.upper() == 'PHILIPS':
            ending = self.type + '_p.dat'
        elif self.device.upper() == 'CARDIOSCREEN':
            ending = self.type + '_c.dat'
        else:
            ValueError('Unknown device name')

        if file_name is None:
            file_loc = Waveform.root + '/Studies/MOLLIE/processed_signals/' + self.volunteer_id + self.study_id + '_pickle/' + ending
        else:
            file_loc = Waveform.root + '/Studies/MOLLIE/processed_signals/' + file_name + '_pickle/' + ending

        saving_var = {}
        for list_item in list:
            exec('saving_var["'+list_item+'"] = self.'+list_item)

        with open(file_loc,'wb') as outf:
            cPickle.dump(saving_var, outf, protocol=pickle.HIGHEST_PROTOCOL)


    def __repr__(self):
        return "General class for waveforms in the MOLLIE dataset"


    def index_t(self, index):
        t_out = self.t[index]
        return t_out

    def index_ts(self, index):
        ts_out = self.ts[index]
        return ts_out

    def get_derivative(self, order=1):
        deriv = self.ts
        for order_num in range(order):
            diff = deriv.diff()
            deriv = diff * self.fs
        return deriv


    def switch_beat_and_ts(self):
        t_holder = self.t
        ts_holder = self.ts
        sqi_holder = self.sqi

        self.t = self.t_beat
        self.ts = self.ts_beat
        self.sqi = self.sqi_beat

        self.t_beat = t_holder
        self.ts_beat = ts_holder
        self.sqi_beat = sqi_holder


    def segment_data(self, t_start, t_end, do_beat_flag = False, verbose_flag = True):

        if do_beat_flag:
            self.switch_beat_and_ts()

        Out_seg = copy.copy(self)
        #Adjust the t_start and t_end

        t_start = max(t_start, self.t[0])
        t_end = min(t_end, self.t[-1])

        if t_start > t_end:
            Out_seg.t = np.array([])
            Out_seg.ts = np.array([])
            Out_seg.sqi = np.array([])
            if do_beat_flag:
                self.switch_beat_and_ts()
                Out_seg.switch_beat_and_ts()
            return Out_seg


        #This function segments the waveform time series at the input times
        index_ts = np.logical_and(self.t >= t_start, self.t <= t_end)
        if not any(index_ts):
            Out_seg.t = np.array([])
            Out_seg.ts = np.array([])
            Out_seg.sqi = np.array([])
            if do_beat_flag:
                self.switch_beat_and_ts()
                Out_seg.switch_beat_and_ts()
            if verbose_flag:
                print('t_start and t_end are not in the signal time vector')
            return Out_seg

        Out_seg.t = Out_seg.t[index_ts]
        Out_seg.t = Out_seg.t - Out_seg.t[0]
        Out_seg.ts = Out_seg.ts[index_ts]
        Out_seg.sqi = Out_seg.sqi[index_ts]

        index_ts = np.where(index_ts)
        if hasattr(self, 'peaks'):
            index_peak = np.where(np.logical_and(self.peaks >= index_ts[0][0], self.peaks <= index_ts[0][-1]))
            Out_seg.peaks = Out_seg.peaks[index_peak[0]]
            Out_seg.peaks -= index_ts[0][0]

        if hasattr(self, 'onsets'):
            index_onsets = np.where(np.logical_and(self.onsets >= index_ts[0][0], self.onsets <= index_ts[0][-1]))
            Out_seg.onsets = Out_seg.onsets[index_onsets[0]]
            Out_seg.onsets -= index_ts[0][0]

        if do_beat_flag:
            self.switch_beat_and_ts()
            Out_seg.switch_beat_and_ts()

        return Out_seg

    def get_good_values(self, sqi_threshold = 0.8, beat_flag =0, nan_flag = False, return_good_idx = False):
        if not beat_flag:
            good = self.ts
            idx = np.where(self.sqi >= sqi_threshold)
            if any(self.sqi < sqi_threshold):
                if nan_flag:
                    good[np.where(self.sqi < sqi_threshold)] = np.nan
                else:
                    good = good[self.sqi > sqi_threshold]
        else:
            good = self.ts_beat
            idx = np.where(self.sqi_beat >= sqi_threshold)
            if any(self.sqi_beat < sqi_threshold):
                if nan_flag:
                    good[np.where(self.sqi_beat < sqi_threshold)] = np.nan
                else:
                    good = good[np.where(self.sqi_beat > sqi_threshold)]

        if return_good_idx:
            return good, idx
        else:
            return good

    def remove_bad_quality_data(self,  sqi_threshold = 0.8, beat_flag =0):

        if not beat_flag:
            self.ts, idx = self.get_good_values(sqi_threshold=sqi_threshold, beat_flag = beat_flag, return_good_idx=True)
            self.t = self.t[idx]
            self.sqi = self.sqi[idx]
        else:
            self.ts_beat, idx = self.get_good_values(sqi_threshold=sqi_threshold, beat_flag = beat_flag, return_good_idx=True)
            self.t_beat = self.t_beat[idx]
            self.sqi_beat = self.sqi_beat[idx]

    def beat_sqi_to_time(self):

        pass

    def get_fft(self, ts= None, fs=None):
        #Should maybe make the time series used a variable -- do we want ts or ts_kalman etc etc
        if ts is None:
            ts = self.ts
        if fs is None:
            fs =self.fs

        #Get length of FFT time series as the next power of two above
        Nfft = 1<<(len(ts)-1).bit_length()

        yf = scipy.fftpack.fft(ts, Nfft)
        N = (Nfft/2) +1
        xf = np.arange(0, N-1) * fs/Nfft
        pw_mag = abs(yf[:int(N-1)] / Nfft)
        sq_mag = pw_mag ** 2
        pw = sq_mag / np.linalg.norm(sq_mag)

        return pw, xf

    def plot_fft(self, annotate= False):
        # function to plot the FFT of a signal

        pw, freq = self.get_fft()


        fig = plt.figure()
        plt.plot(freq, pw, zorder=1)
        plt.grid()

        if annotate:
            f_max = freq[np.argmax(pw)]
            HR_max = 60 * f_max
            print(max(pw))
            label = f'Max freq = {f_max :.2f} Hz. ({HR_max :.2f} bpm)'
            plt.scatter(f_max, max(pw), c='r', zorder=2, label=label)
            plt.legend()




    def get_signal_at_cuff_times(self, BP, align_flag, invert_flag = 0, max_lag = 5*60,
                                 window_size = 20, sqi_threshold = 0.8,
                                 good_beat_threshold = 0.5, debug = False, debug_deep = False,
                                 return_correlation_flag = False, return_delay_flag = False):
        delay = 0
        if align_flag:
            if not hasattr(BP, 'align'):
                BP.get_BP_align()
            BP_align = BP.align

            #Need to align the BP signal with the signal
            delay  = find_delay_signals(Signal_one=BP_align, Signal_two=self,invert_two_flag = invert_flag, max_lag=max_lag, debug =debug_deep)
            delay = delay/ self.fs
            self.t = self.t + delay

        #Now that the signals are aligned we need to get signal values at the time of cuff inflation
        data = []
        num_inflations = len(BP.ts)

        for inflation_idx in range(num_inflations):
            if BP.sqi[inflation_idx] < sqi_threshold:
                data.append([np.nan, np.nan])
                continue
            t_start = BP.t[inflation_idx] - window_size/2
            t_end = BP.t[inflation_idx] + window_size / 2
            #if start and end are not in then skip
            if t_end < self.t[0] or t_start > self.t[-1]:
                data.append([BP.ts[inflation_idx], np.nan])
                # data.append([BP.ts[inflation_idx], 'Error'])
                continue

            signal_seg = self.segment_data(t_start=t_start, t_end=t_end)
            #Check if more than good_beat_threshold proportion are of good quality
            if len(signal_seg.get_good_values())/len(signal_seg.ts) < good_beat_threshold:
                data.append([BP.ts[inflation_idx], np.nan])
                continue

            data.append([BP.ts[inflation_idx], np.nanmean(signal_seg.get_good_values())])


        data = np.array(data)
        if return_correlation_flag:
            is_not_nan = ~(np.isnan(data).any(axis=1))
            correlation = np.corrcoef(data[is_not_nan,0],data[is_not_nan,1] )[0,1]

        if debug:
            sns.regplot(y=data[is_not_nan,0], x=data[is_not_nan,1])
            plt.grid()
            plt.title(f'Correlation coefficient = {correlation:.2f} ')
            plt.xlabel(self.type + ' (' + self.units + ')')
            plt.ylabel('BP (mmHg)')




        if return_correlation_flag and return_delay_flag:
            return data, correlation, delay
        elif return_delay_flag and ~return_correlation_flag:
            return data, delay
        elif ~return_delay_flag and return_correlation_flag:
            return data, correlation
        else:
            return data

    def copy(self):
        return copy.deepcopy(self)



    def show(self, mins_flag = 0):
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        if not mins_flag:
            ax1.plot(self.t, self.ts)
            ax1.set_xlabel('Time (seconds)')
            ax1.set_ylabel(self.type + '(' + self.units + ')', fontsize=12)

            ax2.plot(self.t, self.sqi)
            ax2.set_xlabel('Time (seconds)')
            ax2.set_ylabel('SQI')
        else:
            ax1.plot(self.t/60, self.ts)
            ax1.set_xlabel('Time (mins)')
            ax1.set_ylabel(self.type + '(' + self.units + ')', fontsize=12)

            ax2.plot(self.t/60, self.sqi)
            ax2.set_xlabel('Time (mins)')
            ax2.set_ylabel('SQI')

    def show_beat_and_ts(self, mins_flag= True):

        factor = 60 if mins_flag else 1

        sns.set_style("darkgrid")
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax1.plot(self.t_beat / factor, self.ts_beat, zorder=1)
        ax1.plot(self.t / factor, self.ts, zorder=2)
        plt.ylabel(self.type + '(' + self.units + ')', fontsize=12)

        ax2.plot(self.t_beat / factor, self.sqi_beat, zorder=1)
        ax2.plot(self.t / factor, self.sqi, zorder=2)

        xlabel_str = 'Time (mins)' if mins_flag else 'Time (s)'

        plt.xlabel(xlabel_str, fontsize=12)
        plt.ylabel('SQI', fontsize=12)

    def show_beat(self, mins_flag= True):

        factor = 60 if mins_flag else 1

        sns.set_style("darkgrid")
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax1.plot(self.t_beat / factor, self.ts_beat, zorder=1)
        plt.ylabel(self.type + '(' + self.units + ')', fontsize=12)

        ax2.plot(self.t_beat / factor, self.sqi_beat, zorder=1)

        xlabel_str = 'Time (mins)' if mins_flag else 'Time (s)'

        plt.xlabel(xlabel_str, fontsize=12)
        plt.ylabel('SQI', fontsize=12)
