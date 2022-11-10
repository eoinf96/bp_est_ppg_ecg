from classes import Waveform

#Allow for loading
#And allow for complete recomputing via ECG or PPG if neccesary
class HR_cls(Waveform):
    def __init__(self,  volunteer_id, device_name='Stowood', study_id = 'a'):
        super().__init__(volunteer_id, load_mat=0 ,  study_id = study_id)
        self.type = 'HR'
        self.units = 'bpm'
        self.device = device_name
        self.study_id = study_id

    def load_mat(self):
        HR = super().load_mat()
        self._load(HR=HR)

    def load_pickle(self, file_name = None):
        HR = super().load_pickle(file_name=file_name)
        self._load(HR=HR)

    def _load(self, HR):
        self.ts = HR['ts']
        self.t = HR['t']
        self.sqi = HR['sqi']
        self.fs = HR['fs']

    def save_pickle(self,list = None, file_name = None):
        if list is None:
            list = ['ts', 't', 'sqi', 'fs']
        super().save_pickle(list=list, file_name=file_name)
