'''
 Released under the GNU General Public License
 Copyright (C) 2021  Eoin Finnegan
 eoin.finnegan@eng.ox.ac.uk

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import numpy as np
import pandas as pd
import copy

class MOLLIE_session:
    """
    A class representing a MOLLIE session, containing demographic data
    for multiple participants and methods for selecting and manipulating
    subsets of participants.Note: the demographic data in this implementation is not the actual
    data for each participant due to data protection.
    """

    def __init__(self):
        """
        Initializes a new MOLLIE_session object by reading in demographic
        data from a CSV file and computing BMI values.
        """
        self.df = pd.read_csv('../data/demographics.csv')
        self.df['BMI'] = self.df['Weight'] / (self.df['Height'] / 100) ** 2

        # Hidden variables
        self._df_base = self.df.copy()

    def select_participants(self, volunteer_list):
        """
        Selects a subset of participants based on a list of volunteer IDs.
        """
        self.df = self.df.iloc[volunteer_list]

    def reset_df(self):
        """
        Resets the participant list to its original state.
        """
        self.df = self._df_base

    @classmethod
    def get_demographics_names(cls):
        """
        Returns a list of column names for the demographic data.
        """
        return ['Age', 'Height', 'BMI', 'Weight', 'Sex']

    def get_copy_without_volunteer(self, v_num):
        """
        Returns a deep copy of the current MOLLIE_session object with
        a specified volunteer removed from the participant list.
        """
        MOLLIE_copy = copy.deepcopy(self)

        # Remove the volunteer from the list
        volunteer_list = np.setdiff1d(self.df['ID'], v_num)
        MOLLIE_copy.select_participants(volunteer_list=volunteer_list)

        return MOLLIE_copy
