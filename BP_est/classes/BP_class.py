from classes import Waveform
from loading_functions import *
import seaborn as sns
import numpy as np

############################# Class for BP
# Consider making this a class for SBP, MAP and DBP -- ts_SBP etc... -- That way it saves on reloading
class BP_cls(Waveform):
    p_global_values = {
        'SBP_a': 8554.672535565684,
        'MAP_a': 5790.443980602483,
        'DBP_a': 1795.714494371641,
        'SBP_b': 255.1406520031287,
        'MAP_b': 1215.474250076286,
        'DBP_b': 376.9390975388363
    }

    def __init__(self,  volunteer_id, BP_string = 'All', study_id = 'a'):
        super().__init__(volunteer_id, load_mat=0, study_id=study_id)
        self.units = 'mmHg'
        self.type = BP_string
        self.study_id = study_id
        self.device = 'Philips'
        if BP_string == 'All':
            self.ts_SBP = []
            self.ts_MAP = []
            self.ts_DBP = []

    def load_mat(self):
        BP = super().load_mat()
        self._load(BP=BP)

    def load_pickle(self, file_name = None):
        BP = super().load_pickle(file_name=file_name)
        self._load(BP=BP)

    def _load(self, BP):
        if self.type == 'All':
            self.ts_SBP = BP['ts_SBP']
            self.ts_MAP = BP['ts_MAP']
            self.ts_DBP = BP['ts_DBP']

            self.t = BP['t']
            self.sqi = np.ones((len(self.ts_SBP),))
            self.fs = BP['fs']
            self.t_new = BP['t_new']
        else:
            self.ts = BP['ts']
            self.t = BP['t']
            self.sqi = np.ones((len(self.ts),))
            self.fs = BP['fs']
            self.t_new = BP['t_new']
            self.name = BP['name']

        if self.record_name == '014a':
            #This data is wrong and is causing bad performance
            self.sqi[6] = 0

    def save_pickle(self, list = None, file_name = None):
        if list is None:
            if self.type == 'All':
                list = ['ts_SBP','ts_MAP', 'ts_DBP', 't', 'sqi', 'fs', 't_new']
            else:
                list = ['ts', 't', 'sqi', 'fs', 't_new', 'name']
        super().save_pickle(list = list, file_name=file_name)



    def get_BP_align(self, fs=60):
        BP_align =self.copy()
        BP_align.compute_smoothing_splines(fs=fs)

        self.align = BP_align

    def set_curr_BP(self, BP_string):
        #Method to simply set one of SBP, MAP or DBP as the current BP
        exec('self.ts = self.ts_' + BP_string)
        self.type = BP_string

    def compute_smoothing_splines_all(self, xx=None, fs=None, p_value=None,
                                  debug_mode=False):
        for bp_string in ['SBP', 'MAP', 'DBP']:
            self.set_curr_BP(bp_string)
            self.compute_smoothing_splines(xx, fs, p_value, debug_mode)
            exec('self.ts_'+bp_string +' = self.ts')
            if bp_string != 'DBP':
                #If we are not on our last iteration then we need to put sqi and t back
                self.t = self.t_orig
                self.sqi = self.sqi_orig




    def compute_smoothing_splines(self, xx=None, fs=None, p_value=None,
                                  debug_mode=False):
        if p_value is None:
            #Then we use the global p_value
            p_value = self.p_global_values[self.type + '_' + self.study_id]
        self.p_value = p_value


        if xx is None:
            if fs is None:
                print('You need to define either xx or fs')
            xx = np.arange(self.t[0], self.t[-1] + 1/fs, 1/fs)
            self.fs = fs


        signal = self.copy()
        signal.remove_bad_quality_data()

        #Function to perform smoothing splines for a cuff blood pressure signals
        y = signal.ts
        t = (signal.t - min(signal.t))/(max(signal.t) - min(signal.t))
        tt = (xx - min(signal.t))/(max(signal.t) - min(signal.t)) # Put tt in the spectrum of t -- any t outside of this will be linearly interpolated
        h = np.diff(t)
        N = len(t)
        # Calculate differences
        Q = np.zeros((N, N-2))

        ind = np.arange(N-2)
        Q[ind, ind] =  1. / h[:- 1]
        Q[ind+1, ind] = -1. / h[:- 1]-1. / h[1:]
        Q[ind+2, ind] = 1. / h[1:]

        T = (np.diag(h[1:-1], -1)+np.diag(2 * (h[0:-1] + h[1:]), 0)+np.diag(h[1:- 1], 1)) / 3

        #If you want to add some weights in then you can do it here -- should look into this
        # B = np.transpose(Q) @ Q + p_value*T
        B = np.matmul(np.transpose(Q), Q) + p_value * T

        # c = np.linalg.lstsq(B, p_value * np.transpose(Q) @ y, rcond=None)[0]
        c = np.linalg.lstsq(B, p_value * np.matmul(np.transpose(Q), y), rcond=None)[0]
        # a = y - (Q @c) / p_value
        a = y - (np.matmul(Q ,c)) / p_value
        c = np.pad(c, (1, 1), 'constant', constant_values=(0, 0))
        d = (c[1:]- c[:-1])/(3 * h)
        d = np.append(d, d[-1])
        b = (a[1:]- a[:-1])/h  - c[:-1]* h - d[:-1]*(h**2)
        b = np.append(b, b[-1])
        #Vertically concatenate
        coeffs = np.c_[d, c, b, a]

        # Now compute the splines
        # y_out = np.zeros((len(tt),1))
        y_out = np.zeros((len(tt), ))
        for jj in range(len(tt)):
            index_along_original_time_series = sum(tt[jj] >= i for i in t) -1
            # Check for negative values or values greater than t_end
            if tt[jj] <=0:
                y_out[jj] = np.polyval(coeffs[0, -2:], tt[jj])
            elif tt[jj] >= 1:
                y_out[jj] = np.polyval(coeffs[-1, -2:], tt[jj] - t[index_along_original_time_series])
            elif tt[jj] == t[index_along_original_time_series]:
                y_out[jj] = coeffs[index_along_original_time_series,3]
            else:
                y_out[jj] = np.polyval(coeffs[index_along_original_time_series, :], tt[jj]-t[index_along_original_time_series])


        self.ts_orig = self.ts
        self.t_orig = self.t
        self.sqi_orig = self.sqi
        self.ts = y_out
        self.t = np.transpose(xx)
        self.sqi = np.ones((len(self.ts),))


        # Set sqi = 0 during periods where there is a big gap for each signal
        start_times = {
            '016a' : 200, '002a': 800, '003a' : 750
        }
        end_times = {
            '016a': 800, '002a': 1200, '003a': 1250
        }
        if self.record_name in start_times:
            index = np.where(np.logical_and(self.t >= start_times[self.record_name], self.t <= end_times[self.record_name]))
            self.sqi[index] = 0


        if debug_mode:
            self.plot_spline()


    def plot_spline(self, mins_flag= True):
        factor = 60 if mins_flag else 1
        xlabel_str =  "Time (mins)" if mins_flag else "Time (secs)"

        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

        sns.set_style("darkgrid")
        ax1.plot(self.t_orig/factor, self.ts_orig,linestyle='--', zorder=1, label='Original cuff data')
        ax1.scatter(self.t_orig/factor, self.ts_orig, s = 15, zorder=2)
        good = self.ts
        good[np.where(self.sqi < 0.8)] = np.nan
        ax1.plot(self.t/factor, good, zorder=3, label='Processed cuff data')
        ax1.set_xlabel(xlabel_str, fontsize=12)
        ax1.set_ylabel(self.type + '(mmHg)', fontsize=12)
        ax1.legend()


        ax2.plot(self.t/factor, self.sqi)
        ax2.set_xlabel(xlabel_str, fontsize=12)
        ax2.set_ylabel('SQI', fontsize=12)



    # def show(self):
    #     #A showing function for blood pressure
