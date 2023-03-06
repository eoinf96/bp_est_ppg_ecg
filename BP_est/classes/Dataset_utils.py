'''
Class to handle regression
--
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

import pandas as pd
import numpy as np
from sklearn.model_selection import LeaveOneGroupOut, StratifiedKFold
from sklearn.preprocessing import StandardScaler

###################################################################
### class to perform BP estimation regression
###################################################################
class datasetHandler():
    def __init__(self, df, feature_names, target_name):
        self.df = df
        self.feature_names = feature_names
        self.target_name = target_name

    def return_X_Y(self, to_numpy=False, relevant_IDs=None):
        if relevant_IDs is None:
            relevant_IDs = self.df['ID'].unique()
        if to_numpy:
            return self.df.loc[self.df['ID'].isin(list(relevant_IDs))][self.feature_names].to_numpy(), self.df.loc[self.df['ID'].isin(list(relevant_IDs))][self.target_name].to_numpy()
        else:
            return self.df.loc[self.df['ID'].isin(list(relevant_IDs))][self.feature_names], self.df.loc[self.df['ID'].isin(list(relevant_IDs))][self.target_name]

    def drop_collinear_features(self, VIF_max = 10):
        (X, _) = self.return_X_Y(to_numpy=True)
        dropping = True

        # Center X
        X = StandardScaler().fit_transform(X)
        while dropping:
            R0 = np.corrcoef(X.T)
            VIF = np.diag(np.linalg.inv(R0))
            if max(VIF) > VIF_max:
                ii = np.argmax(VIF)
                X = np.delete(X, ii, 1)
                self.feature_names.remove(self.feature_names[ii])
            else:
                dropping = False



class BP_est_dataset():
    '''
      Arguments
          dataset -- dataset class containing information on the train and test population
          df -- dataframe of feature values
          model_name -- name of the regression model: Linear, RF and SVM supported
          eval_method -- Set to 'CV' to perform Leave One Subject Out Cross Validation. Else Train Test split is used
          add_demographics -- Flag on whether to add demographics to the dataset. 1 for yes. 0 for no
      '''

    # pylint: disable=too-many-instance-attributes
    # pylint: disable=too-many-arguments
    def __init__(
        self,
        study=None,
        df=None,
        df_aug=None,
        BP_name = 'SBP',
        add_demographics=True,
    ):

        self.BP_name = BP_name
        self.demo_added = add_demographics

        # Add df and dataset
        self.Original = None
        self.Augmented = None
        self.waveform_feature_names = None
        if study is not None:
            self.upload_study(study=study)
            if df is not None:
                self.upload_df(df=df, df_aug=df_aug)

        self.ID_i = None

    ##############
    # DF Initialisation
    ##############
    def upload_df(self, df, df_aug):
        '''Upload dataframe containing features here.'''
        df_attr_name = {0: 'Original', 1: 'Augmented'}
        for (ii, curr_df) in enumerate([df, df_aug]):
            if self.demo_added:
                # df_demo = self.dataset.get_demographics_df()
                df_demo = self.study.df.copy()
                df_demo['Sex'] = df_demo['Sex'].map({'F':1, 'M':0})
                # Inner join the df with df_demo
                curr_df = pd.merge(curr_df, df_demo, on = 'ID').dropna()
            else:
                curr_df = curr_df.dropna()

            feature_names =[x for x in list(curr_df.columns) if x not in ['BP', 'SBP', 'MAP', 'DBP', 'ID', 'calibration']]
            curr_df = datasetHandler(df=curr_df,feature_names=feature_names, target_name=self.BP_name)
            setattr(self, df_attr_name[ii], curr_df)
        self.feature_names = self.Original.feature_names
        self.waveform_feature_names = [x for x in feature_names if x not in self.study.get_demographics_names()]

    def upload_study(self, study):
        '''Upload dataset containing demographics.'''
        self.study = study

    def calibrate_dataset(self, std_flag=False, num_cali=5):
        for (ii, curr_df) in enumerate([self.Original.df, self.Augmented.df]):
            '''Subtract norm per feature per subject'''
            if self.demo_added:
                demo_names = self.study.get_demographics_names()
                # Drop demo for now and then re apply them
                df_demo = curr_df.reindex(columns=['ID'] + demo_names)

                # not_demo = np.setdiff1d(list(df.columns), demo_names)
                not_demo = [feat for feat in list(curr_df.columns) if feat not in demo_names]
                curr_df = curr_df.reindex(columns=not_demo)

            unique_ids = curr_df['ID'].unique()
            if 'calibration' not in  list(curr_df.columns):
                if num_cali is not None:
                    curr_df['calibration'] = np.zeros((len(curr_df),))
                    for id in unique_ids:
                        curr_df['calibration'][np.where(curr_df['ID'] == id)[0][:num_cali]] = 1
                else:
                    ValueError('Unable to get the calibration indices. Pass a calibration column or num_cali')
            means = curr_df.loc[curr_df['calibration'] == 1].groupby("ID", as_index=False).mean()
            means = means.reindex(means.index.repeat(curr_df.groupby('ID')['SBP'].count())).reset_index(drop=True)
            stds = curr_df.loc[curr_df['calibration'] == 1].groupby("ID").std()
            stds = stds.reindex(stds.index.repeat(curr_df.groupby('ID')['SBP'].count()))

            means = means.replace(0, np.finfo(float).eps)
            stds = stds.replace(0, np.finfo(float).eps)
            if std_flag:
                curr_df[self.waveform_feature_names] = (curr_df[self.waveform_feature_names] - means[self.waveform_feature_names])/stds[self.waveform_feature_names]
                curr_df[self.BP_name] = (curr_df[self.BP_name] - means[self.BP_name]) / stds[self.BP_name]
            else:
                curr_df[self.waveform_feature_names] = (curr_df[self.waveform_feature_names] - means[self.waveform_feature_names]) / means[self.waveform_feature_names]
                curr_df[self.BP_name] = curr_df[self.BP_name] - means[self.BP_name]

            if self.demo_added: curr_df = pd.merge(curr_df, df_demo, on='ID')
            if ii == 0:
                self.Original.df = curr_df
            else:
                self.Augmented.df = curr_df

    def drop_collinear_features(self, VIF_max=10):
        self.Original.drop_collinear_features(VIF_max=VIF_max)
        self.feature_names = self.Original.feature_names
        self.Augmented.feature_names = self.Original.feature_names

    def ID_iterator(self):
        # Assert no ID!!
        unique_ids = self.df['ID'].unique()
        if self.ID_i is None:
            self.ID_i =unique_ids[0]
        else:
            ii =np.where(unique_ids == self.ID_i)[0][0]
            self.ID_i = unique_ids[ii+1]

    def return_LOSO_train_test_data(self):
        pass


    def generate_cv(self, eval_method='LOSO', nFolds = 5):
        '''Initialise CV Generator using the dataset'''
        (X, y) = self.return_X_Y()
        if eval_method == 'LOSO':
            self.cv = LeaveOneGroupOut().split(
                X=X,
                y=y,
                groups=self.df['ID'],
            )
        elif eval_method == 'CV':
            self.cv = StratifiedKFold(n_splits=nFolds).split(
                X=X,
                y=y,
                groups=self.df['ID'],
            )