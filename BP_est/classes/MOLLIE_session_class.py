import numpy as np
import pandas as pd
# from set_root import set_root
import copy


class MOLLIE_session():

    # Define constants
    __const_volunteer_list = np.arange(1, 32)
    __const_Age = np.array([30,32,31,33,54,32,25,24,27,28,37,47,24,25,30,21,18,28,35,58,26,30,35,28,22,25,43,21,33,23,34])
    __const_Weight = np.array([68, 83, 54, 81, 72, 83, 53, 63, 56, 67, 51, 73, 75, 82, 52, 54, 39, 72, 85, 57, 74, 54, 62, 60,
                        79, 65, 55, 59, 66, 73, 80])
    __const_Height = np.array([174, 180, 161, 173, 168, 184, 162, 159, 167, 187, 157, 161, 173, 185, 162, 162, 152, 178, 180,
                        166, 175, 160, 167, 169, 171, 175, 171, 173, 181, 186, 162])  # Add heights in here
    __const_BMI = np.array([22.5, 25.4,20.83,27.1,25.5,24.5,20.2,24.9,20.1,19.2,20.5,28.2,25.1,23.9,19.9,20.6,16.9,22.7,26.2,20.5,24.2,21.1,22.2,21,26.8,21.2,18.8,19.6,20.1,21.1,30.5])
    __const_Male = np.array([1,1,0,1,1,1,0,0,1,1,0,0,1,1,0,0,0,1,1,0,1,0,0,0,1,1,0,0,1,1,0])
    #These should possibly include 1, 20 and 26
    __const_TrainID = np.array([4, 8, 9, 10, 11, 12, 13, 18, 24, 25, 27, 28, 29, 31])
    __const_TestID = np.array([2, 3, 5, 6, 7, 14, 15, 16, 17, 19, 21, 22, 23, 30])

    def __init__(self, study_id = 'a'):
        self.study_id = study_id
        self.reset_volunteer_list()


    def set_volunteer_list(self, volunteer_list, remove_noisy_volunteers_flag = 1):
        self.volunteer_list = np.array(volunteer_list)
        if remove_noisy_volunteers_flag:
            # if self.study_id == 'a':
            #     remove = np.array([1, 20, 26, 27])
            # elif self.study_id =='b':
            #     remove = np.array([1, 2, 9])
            # else:
            #     raise ValueError('Unknown study ID')
            remove = np.array([1, 2, 3, 9, 20, 26])
            self.volunteer_list = np.setdiff1d(self.volunteer_list, remove)

        self.num_volunteers = len(self.volunteer_list)
        self.Age = self.__const_Age[self.volunteer_list -1]
        self.Height = self.__const_Height[self.volunteer_list -1]
        self.Weight = self.__const_Weight[self.volunteer_list -1]
        self.BMI = self.__const_BMI[self.volunteer_list-1]
        self.Male = self.__const_Male[self.volunteer_list-1]
        # self.Fitzpatrick
        self.TrainID = np.intersect1d(self.__const_TrainID, self.volunteer_list)
        self.TestID = np.intersect1d(self.__const_TestID, self.volunteer_list)
        self.root = set_root()
        self.figure_path = self.root + '/figures/MOLLIE/'
        self.results_path = self.root + '/Studies/MOLLIE/results/'

    def reset_volunteer_list(self):
        self.set_volunteer_list(volunteer_list=self.__const_volunteer_list)

    def get_volunteer_id_name(self, volunteer_num):
        if isinstance(volunteer_num, np.ndarray):
            volunteer_id = [ '00' + str(s)  for s in volunteer_num]
            volunteer_id = [s[-3:] for s in volunteer_id]
        else:
            volunteer_id = '00' + str(volunteer_num)
            volunteer_id = volunteer_id[-3:]

        return volunteer_id


    def get_volunteer_id_name_loop_idx(self, v_idx):
        volunteer_id = self.volunteer_list[v_idx-1]
        volunteer_id = '00' + str(volunteer_id)
        volunteer_id = volunteer_id[-3:]

        return volunteer_id


    def get_demographics_df(self):

        a = {
            'Age': self.Age,
            'BMI': self.BMI,
            'Weight': self.Weight,
            'Gender' : self.Male,
            'ID': self.get_volunteer_id_name(self.volunteer_list)
        }
        df_demo = pd.DataFrame(a, columns=['Age', 'BMI', 'Weight','Gender', 'ID'])
        df_demo.index = df_demo['ID']
        df_demo = df_demo.drop(columns='ID')

        return df_demo

    @classmethod
    def get_demographics_names(cls):
        return ['Age', 'BMI', 'Weight', 'Gender', 'ID']


    def get_copy_without_volunteer(self, v_num):
        MOLLIE_copy = copy.deepcopy(self)

        #Remove the v_nums from the list
        volunteer_list = np.setdiff1d(self.volunteer_list, v_num)
        MOLLIE_copy.set_volunteer_list(volunteer_list=volunteer_list)

        return MOLLIE_copy
