import numpy as np
import pandas as pd
# from set_root import set_root
import copy


class MOLLIE_session():

    # NOTE THAT DUE TO DATA PROTECTION THESE DEMOGRAPHIC CONSTANTS HAVE BEEN CHANGED AND ARE NOT!! THE ACTUAL DEMOGRAPHICS FOR EACH PARTICIPANT

    def __init__(self):
        self.df = pd.read_csv('data/demographics.csv')
        self.df['BMI'] = self.df['Weight'] / (self.df['Height'] / 100) ** 2

        # Hidden variables
        self._df_base = self.df.copy()

    def select_participants(self, volunteer_list):
        self.df = self.df.iloc[volunteer_list]

    def reset_df(self):
        self.df = self._df_base

    @classmethod
    def get_demographics_names(cls):
        return ['Age', 'Height', 'BMI', 'Weight', 'Sex']


    def get_copy_without_volunteer(self, v_num):
        MOLLIE_copy = copy.deepcopy(self)

        #Remove the v_nums from the list
        volunteer_list = np.setdiff1d(self.df['ID'], v_num)
        MOLLIE_copy.select_participants(volunteer_list=volunteer_list)

        return MOLLIE_copy
